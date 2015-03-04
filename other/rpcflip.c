// RPCFLIP
//
// A "Nadir" viewer for pleiades images, based on vflip by Gabriele Facciolo
//
// c99 -O3 rpcflip.c iio.o -o rpcflip -lX11 -ltiff -lm
//
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "iio.h"

#define TIFFU_OMIT_MAIN
#include "tiffu.c"

#define DONT_USE_TEST_MAIN
#include "rpc.c"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.c"

#include "xmalloc.c"

// TODO:
//
// 1. mouse scroll (with fast panning skipping the pixel-wise rpc)
// 2. transparent detection of the 12 significant bits on the 16 bit smaples
// 3. automatic contrast adjust à la qauto
// 4. local interpolation
// 5. setup base_h automatically using the SRTM4
// 6. draw epipolar curves
// 7. load an arbitrary DEM
// 8. view the DEM (and the SRTM4)
//




#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EARTH_RADIUS 6378000.0

#define WHEEL_FACTOR 1.4


// data for a single view
struct pan_view {
	// main image data
	int w, h;                  // pixel dimensions of the P image
	struct tiff_octaves tg[1]; // P
	struct tiff_octaves tc[1]; // MS

	// RGB preview
	uint8_t *preview;
	int pw, ph;

	// calibration
	struct rpc r[1]; 
	double P[8]; // affine approximation to the RPC
	double rgbiox, rgbioy; // offset between P and MS (in P pixel units)
};


#define MAX_VIEWS 30
struct pan_state {
	// 1. data for each view
	int nviews;
	int current_view;
	struct pan_view view[MAX_VIEWS];

	// 2. geographic coordinates (nominal pixel lon/lat resolution)
	double base_h;
	double lon_0, lon_d;
	double lat_0, lat_d;

	// 3. view port parameters
	double a, b; // linear contrast change
	double zoom_factor, offset_x, offset_y; // window <-> geography
	int octave; // (related to zoom-factor)

	// 4. visualization options
	int fix_viewports;
	int force_exact;
};

static uint8_t float_to_uint8(float x)
{
	if (x < 1) return 0;
	if (x > 254) return 255;
	return x;
}

static uint8_t *load_nice_preview(char *fg, char *fc, int *w, int *h)
{
	int wg, hg, wc, hc;
	uint8_t *g = (void*)iio_read_image_uint8_rgb(fg, &wg, &hg);
	uint8_t *c = (void*)iio_read_image_uint8_rgb(fc, &wc, &hc);
	uint8_t *r = xmalloc(3 * wg * hg);
	memset(r, 0, 3 * wg * hg);

	for (int j = 0; j < hg; j++)
	for (int i = 0; i < wg; i++)
	if (i < wc && j < hc)
	{
		int cidx = j*wc + i;
		int gidx = j*wg + i;
		float R = c[3*cidx + 0];
		float G = c[3*cidx + 1];
		float B = c[3*cidx + 2];
		float nc = hypot( R, hypot(G, B));
		float P = g[3*gidx];
		r[3*gidx + 0] = float_to_uint8(3 * P * R / nc);
		r[3*gidx + 1] = float_to_uint8(3 * P * G / nc);
		r[3*gidx + 2] = float_to_uint8(3 * P * B / nc);
	}

	//iio_save_image_uint8_vec("/tmp/prev_g.png", g, wg, hg, 3);
	//iio_save_image_uint8_vec("/tmp/prev_c.png", c, wc, hc, 3);
	//iio_save_image_uint8_vec("/tmp/prev_r.png", r, wg, hg, 3);

	free(g);
	free(c);
	*w = wg;
	*h = hg;
	return r;
}

static void init_view(struct pan_view *v,
		char *fgtif, char *fgpre, char *fgrpc, char *fctif, char *fcpre)
{
	// build preview (use only color, by now)
	//v->preview = (void*)iio_read_image_uint8_rgb(fgpre, &v->pw, &v->ph);
	//fprintf(stderr, "loaded preview image from file \"%s\" (%d %d)\n",
	//		fcpre, v->pw, v->ph);
	v->preview = load_nice_preview(fgpre, fcpre, &v->pw, &v->ph);

	// load P and MS
	int megabytes = 400;
	tiff_octaves_init(v->tg, fgtif, megabytes);
	tiff_octaves_init(v->tc, fctif, megabytes);
	v->w = v->tg->i->w;
	v->h = v->tg->i->h;
	v->rgbiox = v->rgbioy = 4; // normally untouched

	// load P's RPC
	read_rpc_file_xml(v->r, fgrpc);
}

static void setup_nominal_pixels_according_to_first_view(struct pan_state *e)
{
	struct pan_view *v = e->view + 0;

	// center of the image in image coordinates
	double cx = v->w / 2;
	double cy = v->h / 2;

	// center of the image in geographic coordinates
	double center[3], csouth[3];
	eval_rpc(center, v->r, cx, cy    , e->base_h);
	eval_rpc(csouth, v->r, cx, cy + 1, e->base_h);

	// stepsize given by 1 pixel to the south
	// XXX WARNING WRONG TODO FIXME : assumes North-South oriented image (!)
	double lat_step = csouth[1] - center[1];
	double latitude = center[1] * (M_PI/180);
	double lonfactor = cos(latitude);
	double lon_step = -lat_step / lonfactor;

	double lon_step_m = (lon_step*lonfactor) * (M_PI/180) * EARTH_RADIUS;
	double lat_step_m = lat_step * (M_PI/180) * EARTH_RADIUS;
	fprintf(stderr, "projection center (%g %g %g) => (%.8lf %.8lf)\n",
			cx, cy, e->base_h, center[0], center[1]);
	fprintf(stderr, "lon_step = %g (%g meters)\n", lon_step, lon_step_m);
	fprintf(stderr, "lat_step = %g (%g meters)\n", lat_step, lat_step_m);

	// fill-in the fields
	e->lat_d = lat_step;
	e->lon_d = lon_step;
	e->lon_0 = center[0];
	e->lat_0 = center[1];
}

static void init_state(struct pan_state *e,
		char *gi[], char *gp[], char *gr[], char *ci[], char *cp[],
		int n)
{
	// set views
	assert(n < MAX_VIEWS);
	e->nviews = n;
	e->current_view = 0;
	for (int i = 0; i < n; i++)
		init_view(e->view + i, gi[i], gp[i], gr[i], ci[i], cp[i]);

	// set geography
	e->base_h = 0;
	e->lon_0 = e->lat_0 = e->lon_d = e->lat_d = NAN;
	setup_nominal_pixels_according_to_first_view(e);


	// set viewport
	e->a = 1;
	e->b = 0;
	e->zoom_factor = 1;
	e->offset_x = e->offset_y = 0;

	// set options
	e->fix_viewports = 1;
	e->force_exact = 0;
	e->octave = 0;
}

static struct pan_view *obtain_view(struct pan_state *e)
{
	assert(e->current_view >= 0);
	assert(e->current_view < MAX_VIEWS);
	assert(e->current_view < e->nviews);
	return e->view + e->current_view;
}

static int obtain_octave(struct pan_state *e)
{
	// The octave depends on the zoom factor.
	// It is not exactly an octave, but an index to an image to query
	//
	// Possible octaves:
	//
	// 0: get the colors from the preview
	// 1: get the colors from the MS image
	// 2: get the intensities from the P image, and hues from the MS
	// 3: get the rgb from P (resulting in a grayscale image)
#ifdef SELECTABLE_OCTAVE
	assert(e->octave <= 3);
	assert(e->octave >= 0);
	return e->octave;
#else
	if (e->zoom_factor < 0.18) return 0;
	if (e->zoom_factor > 0.5) return 2;
	return 1;
#endif
}

// compute an affine approximation of the RPC function
static void approximate_projection(double P[8], struct rpc *R,
		double x, double y, double h)
{
	// finite difference steps
	double eps_L = 0.000001; // (in geographic DEGREES)
	double eps_H = 1; // (in meters)

	// eval the function at a coordinate tetrahedron
	double v0[3], vx[3], vy[3], vh[3];
	eval_rpci(v0, R, x         , y         , h        );
	eval_rpci(vx, R, x + eps_L , y         , h        );
	eval_rpci(vy, R, x         , y + eps_L , h        );
	eval_rpci(vh, R, x         , y         , h + eps_H);

	// compute forward differences
	double pp = v0[0];
	double qq = v0[1];
	double px = ( vx[0] - v0[0] ) / eps_L;
	double py = ( vy[0] - v0[0] ) / eps_L;
	double ph = ( vh[0] - v0[0] ) / eps_H;
	double qx = ( vx[1] - v0[1] ) / eps_L;
	double qy = ( vy[1] - v0[1] ) / eps_L;
	double qh = ( vh[1] - v0[1] ) / eps_H;

	// fill-in the coefficients of the affine approximation
	P[0]=px; P[1]=py; P[2]=ph;
	P[4]=qx; P[5]=qy; P[6]=qh;
	P[3] = pp - px * x - py * y - ph * h;
	P[7] = qq - qx * x - qy * y - qh * h;
}

static void apply_projection(double out[3], double P[8], double x[3])
{
	out[0] = P[0] * x[0] + P[1] * x[1] + P[2] * x[2] + P[3];
	out[1] = P[4] * x[0] + P[5] * x[1] + P[6] * x[2] + P[7];
	out[2] = x[2];
}

static void window_to_raster(double xy[2], struct pan_state *e, int i, int j)
{
	xy[0] = e->offset_x + i / e->zoom_factor;
	xy[1] = e->offset_y + j / e->zoom_factor;
}

static void window_to_geo(double lonlat[2], struct pan_state *e, int i, int j)
{
	//fprintf(stderr, "\twindow_to_geo(%d %d)\n", i, j);
	double x = e->offset_x + i / e->zoom_factor;
	double y = e->offset_y + j / e->zoom_factor;
	//fprintf(stderr, "\twindow_to_geo xy = (%g %g)\n", x, y);
	lonlat[0] = e->lon_0 + x * e->lon_d;
	lonlat[1] = e->lat_0 + y * e->lat_d;
	//fprintf(stderr, "\twindow_to_geo lonlat = (%g %g)\n", lonlat[0], lonlat[1]);
}

static void window_to_image_ap(double out[2], struct pan_state *e, int i, int j)
{
	double lonlath[3];
	window_to_geo(lonlath, e, i, j);
	lonlath[2] = e->base_h;

	struct pan_view *v = obtain_view(e);
	double p[3];
	apply_projection(p, v->P, lonlath);
	out[0] = p[0];
	out[1] = p[1];
}

static void window_to_image_ex(double out[2], struct pan_state *e, int i, int j)
{
	double lonlath[3];
	window_to_geo(lonlath, e, i, j);
	lonlath[2] = e->base_h;

	struct pan_view *v = obtain_view(e);
	double p[3];
	eval_rpci(p, v->r, lonlath[0], lonlath[1], e->base_h);

	out[0] = p[0];
	out[1] = p[1];
}

static float rgb_getsamplec(void *vx, int w, int h, int i, int j, int l)
{
	int pd = 3;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	unsigned char (*x)[w][pd] = vx;
	return x[j][i][l];
}

static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static void preview_at(float *result, struct pan_view *v, double p, double q)
{
	int ip = p;
	int iq = q;
	int w = v->pw;
	int h = v->ph;
	unsigned char *x = v->preview;
	for (int l = 0; l < 3; l++) {
		float a = rgb_getsamplec(x, w, h, ip  , iq  , l);
		float b = rgb_getsamplec(x, w, h, ip+1, iq  , l);
		float c = rgb_getsamplec(x, w, h, ip  , iq+1, l);
		float d = rgb_getsamplec(x, w, h, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

static int insideP(int w, int h, int x, int y)
{
	return x >= 0 && y >= 0 && x < w && y < h;
}

static double getgreen(double c[4])
{
	return c[1] * 0.6 + c[3] * 0.2;
}

static void pixel_from_preview(float *r, struct pan_view *v, double p, double q)
{
	int ip = p * v->pw / v->w;
	int iq = q * v->ph / v->h;
	float fp = 1 * p * v->pw / v->tg->i->w;
	float fq = 1 * q * v->ph / v->tg->i->h;
	if (insideP(v->pw, v->ph, ip, iq)) {
		preview_at(r, v, fp, fq);
	} else {
		r[0] = 200;
		r[1] = 50;
		r[2] = 100;
	}
}

static void bilinear_ms(float *result, struct pan_view *v, double p, double q)
{
	int ip = p;
	int iq = q;
	int w = v->pw;
	int h = v->ph;
	unsigned char *x = v->preview;
	for (int l = 0; l < 3; l++) {
		float a = rgb_getsamplec(x, w, h, ip  , iq  , l);
		float b = rgb_getsamplec(x, w, h, ip+1, iq  , l);
		float c = rgb_getsamplec(x, w, h, ip  , iq+1, l);
		float d = rgb_getsamplec(x, w, h, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

static void pixel_from_ms(float *out, struct pan_view *v, double p, double q)
{
	int fmt = v->tg->i->fmt;
	int bps = v->tg->i->bps;
	int spp = v->tg->i->spp;
	int ss = bps / 8;
	int ipc = (p + v->rgbiox) / 4; //(p + e->rgbiox) / factor;
	int iqc = (q + v->rgbioy) / 4; //(q + e->rgbioy) / factor;
	char *pix_rgbi = tiff_octaves_getpixel(v->tc, 0, ipc, iqc);

	double c[4];
	c[0] = from_sample_to_double(pix_rgbi + 0*ss, fmt, bps);
	c[1] = from_sample_to_double(pix_rgbi + 1*ss, fmt, bps);
	c[2] = from_sample_to_double(pix_rgbi + 2*ss, fmt, bps);
	c[3] = from_sample_to_double(pix_rgbi + 3*ss, fmt, bps);
	c[1] = getgreen(c);

	double g = (c[0] + c[1] + c[2] + c[3]) / 4;

	double nc = 4*hypot(c[0], hypot(c[1], c[2]));
	out[0] = c[0] * g / nc;
	out[1] = c[1] * g / nc;
	out[2] = c[2] * g / nc;
}

static void pixel_from_pms(float *out, struct pan_view *v, double p, double q)
{
	int fmt = v->tg->i->fmt;
	int bps = v->tg->i->bps;
	int spp = v->tg->i->spp;
	int ss = bps / 8;
	int ipg = p;
	int iqg = q;
	int ipc = (p + v->rgbiox) / 4;
	int iqc = (q + v->rgbioy) / 4;
	//char *pix_gray = tiff_octaves_getpixel_shy(v->tg, 0, ipg, iqg);
	char *pix_gray = tiff_octaves_getpixel(v->tg, 0, ipg, iqg);
	char *pix_rgbi = tiff_octaves_getpixel(v->tc, 0, ipc, iqc);

	double c[4];
	c[0] = from_sample_to_double(pix_rgbi + 0*ss, fmt, bps);
	c[1] = from_sample_to_double(pix_rgbi + 1*ss, fmt, bps);
	c[2] = from_sample_to_double(pix_rgbi + 2*ss, fmt, bps);
	c[3] = from_sample_to_double(pix_rgbi + 3*ss, fmt, bps);
	c[1] = getgreen(c);
	double g;
	if (pix_gray)
		g = from_sample_to_double(pix_gray, fmt, bps);
	else
		g = (c[0] + c[1] + c[2] + c[3]) / 4;
	double nc = 4*hypot(c[0], hypot(c[1], c[2]));
	out[0] = c[0] * g / nc;
	out[1] = c[1] * g / nc;
	out[2] = c[2] * g / nc;
}

static void pixel_from_p(float *out, struct pan_view *v, double p, double q)
{
	pixel_from_p(out, v, p, q);
	int fmt = v->tg->i->fmt;
	int bps = v->tg->i->bps;
	int spp = v->tg->i->spp;
	int ss = bps / 8;
	int ipc = p;
	int iqc = q;
	char *pix_gray = tiff_octaves_getpixel(v->tg, 0, ipc, iqc);

	double g = from_sample_to_double(pix_gray, fmt, bps)/4;
	out[0] = g;
	out[1] = g;
	out[2] = g;
}

// evaluate the value a position (p,q) in image coordinates
static void pixel(float *out, struct pan_view *v, double p, double q, int o)
{
	if (p < 0 || q < 0 || p >= v->w || q >= v->w) {
		out[0] = 0;
		out[1] = 255;
		out[2] = 0;
	}
	if      (o == 0) pixel_from_preview (out, v, p, q);
	else if (o == 1) pixel_from_ms      (out, v, p, q);
	else if (o == 2) pixel_from_pms     (out, v, p, q);
	else if (o == 3) pixel_from_p       (out, v, p, q);
	else exit(fprintf(stderr,"ERROR: bad octave %d\n", o));
}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->zoom_factor;
	e->offset_y -= dy/e->zoom_factor;

	f->changed = 1;
}

static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double c[2];
	window_to_raster(c, e, x, y);

	e->zoom_factor *= F;
	e->offset_x = c[0] - x/e->zoom_factor;
	e->offset_y = c[1] - y/e->zoom_factor;
	fprintf(stderr, "\t zoom changed to %g\n", e->zoom_factor);

	f->changed = 1;
}

static void action_multiply_contrast(struct FTR *f, double fac)
{
	struct pan_state *e = f->userdata;
	e->a *= fac;
	fprintf(stderr, "\t constrast factor changed to %g\n", e->a);
	f->changed = 1;
}

static void action_offset_base_h(struct FTR *f, double d)
{
	struct pan_state *e = f->userdata;

	e->base_h += d;

	fprintf(stderr, "base_h = %g\n", e->base_h);

	f->changed = 1;
}

static void action_increase_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, WHEEL_FACTOR);
}

static void action_decrease_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, 1.0/WHEEL_FACTOR);
}

static void action_toggle_exact_rpc(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->force_exact = !e->force_exact;
	if (e->force_exact)
		fprintf(stderr, "use EXACT calibration\n");
	else
		fprintf(stderr, "use APPROXIMATE calibration\n");
	f->changed = 1;
}

static void action_select_view(struct FTR *f, int i)
{
	struct pan_state *e = f->userdata;
	if (i >= 0 && i < e->nviews)
	{
		fprintf(stderr, "selecting view %d\n", i);
		e->current_view = i;
		f->changed = 1;
	}
}

static void action_select_octave(struct FTR *f, int i)
{
	struct pan_state *e = f->userdata;
	if (i >= 0 && i < 4)
	{
		fprintf(stderr, "selecting octave %d\n", i);
		e->octave = i;
		f->changed = 1;
	}
}

static void action_offset_rgbi(struct FTR *f, double dx, double dy)
{
	struct pan_state *e = f->userdata;
	struct pan_view *v = obtain_view(e);
	v->rgbiox += dx;
	v->rgbioy += dy;
	fprintf(stderr, "rgbi offset chanted to %g %g\n", v->rgbiox, v->rgbioy);
	f->changed = 1;
}


// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		unsigned int pu = p;
		r = pu - (-n) % pu;
		if (r == pu)
			r = 0;
	}
	return r;
}

static void action_cycle_view(struct FTR *f, int d)
{
	struct pan_state *e = f->userdata;
	e->current_view = good_modulus(e->current_view + d, e->nviews);
	f->changed = 1;
}

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;
	struct pan_view *v = obtain_view(e);

	fprintf(stderr, "pan exposer\n");

	// compute 
	double ll[3];
	window_to_geo(ll, e, f->w/2, f->h/2);
	approximate_projection(v->P, v->r, ll[0], ll[1], e->base_h);

	int o = obtain_octave(e);

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		if (e->force_exact)
			window_to_image_ex(p, e, i, j);
		else
			window_to_image_ap(p, e, i, j);
		float c[3];
		pixel(c, v, p[0], p[1], o);
		unsigned char *cc = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < 3; l++)
		{
			float g = e->a * c[l] + e->b;
			if      (g < 0)   cc[l] = 0  ;
			else if (g > 255) cc[l] = 255;
			else              cc[l] = g  ;
		}
	}
	f->changed = 1;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	////fprintf(stderr, "button b=%d m=%d\n", b, m);
	//if (b == FTR_BUTTON_UP && m == FTR_MASK_SHIFT) {
	//	action_contrast_span(f, 1/1.3); return; }
	//if (b == FTR_BUTTON_DOWN && m == FTR_MASK_SHIFT) {
	//	action_contrast_span(f, 1.3); return; }
	//if (b == FTR_BUTTON_RIGHT && m == FTR_MASK_CONTROL) {
	//	action_reset_zoom_only(f, x, y); return; }
	////if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	//if (b == FTR_BUTTON_UP && m & (FTR_MASK_SHIFT|FTR_MASK_CONTROL)) {
	//	action_offset_base_h(f, +10); return; }
	//if (b == FTR_BUTTON_DOWN && m & (FTR_MASK_SHIFT|FTR_MASK_CONTROL)) {
	//	action_offset_base_h(f, -10); return; }

	if (b == FTR_BUTTON_UP && m & FTR_MASK_CONTROL) {
		action_cycle_view(f, +1); return; }
	if (b == FTR_BUTTON_DOWN && m & FTR_MASK_CONTROL) {
		action_cycle_view(f, -1); return; }

	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	//if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	static double ox = 0, oy = 0;

	if (m == FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	//if (m == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	//if (m == FTR_MASK_SHIFT)    action_center_contrast_at_point(f, x, y);
	//if (m == FTR_MASK_CONTROL)  action_base_contrast_at_point(f, x, y);

	ox = x;
	oy = y;
}

//void key_handler_print(struct FTR *f, int k, int m, int x, int y)
//{
//	fprintf(stderr, "key pressed %d '%c' (%d) at %d %d\n",
//			k, isalpha(k)?k:' ', m, x, y);
//}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isalnum(k)?k:' ', m, x, y);

	if (k == '+' || k == '=') action_increase_zoom(f, f->w/2, f->h/2);
	if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);
	//if (k == '+') action_change_zoom_by_factor(f, f->w/2, f->h/2, 2);
	//if (k == '-') action_change_zoom_by_factor(f, f->w/2, f->h/2, 0.5);
	//if (k == 'p') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.1);
	//if (k == 'm') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.1);
	//if (k == 'P') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.006);
	//if (k == 'M') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.006);
	
	if (k == 'j') action_offset_rgbi(f, 0, -1);
	if (k == 'k') action_offset_rgbi(f, 0, 1);
	if (k == 'h') action_offset_rgbi(f, 1, 0);
	if (k == 'l') action_offset_rgbi(f, -1, 0);
	//if (k == ' ') action_cycle(f, 1);
	//if (k == '\b') action_cycle(f, -1);

	//if (isdigit(k)) action_select_image(f, (k-'0')?(k-'0'-1):10);
	if (k == 'a') action_multiply_contrast(f, 1.3);
	if (k == 's') action_multiply_contrast(f, 1/1.3);
	//if (k == 'b') action_contrast_change(f, 1, 1);
	//if (k == 'B') action_contrast_change(f, 1, -1);
	//
	if (k == 'u') action_offset_base_h(f, +10);
	if (k == 'd') action_offset_base_h(f, -10);

	if (k == 'e') action_toggle_exact_rpc(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// image flip operations
	if (k == ' ') action_cycle_view(f, 1);
	if (k == '\b') action_cycle_view(f, -1);

	if (isdigit(k)) action_select_octave(f, (k-'0')?(k-'0'-1):10);

	// arrows move the viewport
	if (k > 1000) {
		int d[2] = {0, 0};
		int inc = -10;
		if (m & FTR_MASK_SHIFT  ) inc /= 10;
		if (m & FTR_MASK_CONTROL) inc *= 10;
		switch (k) {
		case FTR_KEY_LEFT : d[0] -= inc; break;
		case FTR_KEY_RIGHT: d[0] += inc; break;
		case FTR_KEY_UP   : d[1] -= inc; break;
		case FTR_KEY_DOWN : d[1] += inc; break;
		}
		if (k == FTR_KEY_PAGE_UP)   d[1] = +f->h/3;
		if (k == FTR_KEY_PAGE_DOWN) d[1] = -f->h/3;
		action_offset_viewport(f, d[0], d[1]);
	}

	//// if 'k', do weird things
	//if (k == 'k') {
	//	fprintf(stderr, "setting key_handler_print\n");
	//	ftr_set_handler(f, "key", key_handler_print);
	//}
}


#define BAD_MIN(a,b) a<b?a:b

int main_pan(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// process input arguments
	char *filename_preview = pick_option(&c, &v, "p", "");
	if (c < 5 || (c - 1) % 5 != 0) {
		fprintf(stderr, "usage:\n\t"
				"%s [P.TIF P.JPG P.RPC MS.TIF MS.JPG]+\n", *v);
		//                0  1     2     3     4      5
		return c;
	}
	int n = (c - 1)/5;
	fprintf(stderr, "we have %d views\n", n);
	char *filename_gtif[n];
	char *filename_gpre[n];
	char *filename_grpc[n];
	char *filename_ctif[n];
	char *filename_cpre[n];
	for (int i = 0; i < n; i++)
	{
		filename_gtif[i] = v[5*i+1];
		filename_gpre[i] = v[5*i+2];
		filename_grpc[i] = v[5*i+3];
		filename_ctif[i] = v[5*i+4];
		filename_cpre[i] = v[5*i+5];
		fprintf(stderr, "view %d\n", i);
		fprintf(stderr, "\tgtif %s\n", filename_gtif[i]);
		fprintf(stderr, "\tgpre %s\n", filename_gpre[i]);
		fprintf(stderr, "\tgrpc %s\n", filename_grpc[i]);
		fprintf(stderr, "\tctif %s\n", filename_ctif[i]);
		fprintf(stderr, "\tcpre %s\n", filename_cpre[i]);
	}

	// start state
	struct pan_state e[1];
	init_state(e,
			filename_gtif,
			filename_gpre,
			filename_grpc,
			filename_ctif,
			filename_cpre,
			n);

	// open window
	struct FTR f = ftr_new_window(320, 320);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , pan_key_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", pan_exposer);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
