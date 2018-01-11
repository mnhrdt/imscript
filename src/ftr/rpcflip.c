// RPCFLIP(1)                  imscript                   RPCFLIP(1)        {{{1
//
// A "Nadir" viewer for pleiades images, based on vflip by Gabriele Facciolo
//
// cc -O3 rpcflip.c iio.o ftr.o egm96.o -o rpcflip -lX11 -ltiff -lm -ljpeg -lpng
//
// TODO:
//
// 1. draw epipolar curves
// 2. load an arbitrary DEM
// 3. view the DEM (and the SRTM4)
// 4. reorganize coordinate changes to simplify the code
//

// #includes {{{1
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h> // for getpid only

#include "tiff_octaves_rw.c"
#include "srtm4o.c"

#define DONT_USE_TEST_MAIN
#include "rpc2.c"

//double srtm4o(double,double,int);
double egm96(double,double);

#include "ftr.h"

#include "iio.h"
#include "xmalloc.c"


// #defines {{{1

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define EARTH_RADIUS 6378000.0
#define WHEEL_FACTOR 2.0
#define BAD_MIN(a,b) ((a)<(b)?(a):(b))


// pan_view and pan_state {{{1

// data for a single view
struct pan_view {
	// main image data
	int w, h;                  // pixel dimensions of the P image
	struct tiff_octaves tg[1]; // P
	struct tiff_octaves tc[1]; // MS
	int gray_only;

	// RGB preview of the whole image
	uint8_t *preview;
	int pw, ph;
	char *pfg, *pfc; // filenames

	// calibration
	struct rpc r[1];  // PAN rpc
	struct rpc rc[1]; // MSI rpc
	double P[8]; // local affine approximation of "raster to image"
	double rgbiox, rgbioy; // offset between P and MS (in P pixel units)

	// current display of this view
	float *fdisplay;
	uint8_t *display;
	int dw, dh;
	int repaint;
};

#define MAX_VIEWS 60
struct pan_state {
	// 1. data for each view
	int nviews;
	int current_view;
	struct pan_view view[MAX_VIEWS];

	// 2. geographic coordinates (nominal pixel lon/lat resolution)
	double lon_0, lon_d;
	double lat_0, lat_d;
	double base_h;

	// 3. view port parameters
	double a, b; // linear contrast change
	double zoom_factor, offset_x, offset_y; // window <-> geography

	// 4. visualization options
	int force_exact;
	int image_space;
	int diff_mode;
	int show_vertdir, vdx, vdy;
	int interpolation_order;
	int image_rotation_status;
	int qauto;
	int log_scale;
	int show_srtm4, so;
	int srtm4_base;
};

// TODO: integrate this variable into the state
static int msoctaves_instead_of_preview;

// generic utility functions {{{1

static uint8_t float_to_uint8(float x)
{
	if (x < 1) return 0;
	if (x > 254) return 255;
	return x;
}

static int insideP(int w, int h, int x, int y)
{
	return x >= 0 && y >= 0 && x < w && y < h;
}


// small image I/O {{{1
//static
uint8_t *load_nice_preview(char *fg, char *fc, int *w, int *h)
{
	//fprintf(stderr, "load nice preview!\n");
	int wg, hg, wc, hc;
	uint8_t *g = (void*)iio_read_image_uint8_rgb(fg, &wg, &hg);
	uint8_t *c = NULL;
	if (fc)
		c = (void*)iio_read_image_uint8_rgb(fc, &wc, &hc);
	else {
		wc = wg;
		hc = hg;
	}
	uint8_t *r = xmalloc(3 * wg * hg);
	memset(r, 0, 3 * wg * hg);

	for (int j = 0; j < hg; j++)
	for (int i = 0; i < wg; i++)
	if (i < wc && j < hc)
	{
		int cidx = j*wc + i;
		int gidx = j*wg + i;
		float R, G, B;
		if (c) {
			R = c[3*cidx + 0];
			G = c[3*cidx + 1];
			B = c[3*cidx + 2];
		} else R = G = B = 127;
		float nc = hypot( R, hypot(G, B));
		float P = g[3*gidx];
		r[3*gidx + 0] = float_to_uint8(3 * P * R / nc);
		r[3*gidx + 1] = float_to_uint8(3 * P * G / nc);
		r[3*gidx + 2] = float_to_uint8(3 * P * B / nc); }

	free(g);
	if (c) free(c);
	*w = wg;
	*h = hg;
	return r;
}

// state initialization functions {{{1
static void init_view(struct pan_view *v,
		char *fgtif, char *fgpre, char *fgrpc, char *fctif, char *fcpre)
{
	// build preview (use only color, by now)
	//v->preview = load_nice_preview(fgpre, fcpre, &v->pw, &v->ph);
	v->preview = NULL;
	v->pfg = fgpre;
	v->pfc = fcpre;
	v->gray_only = 0;

	// load P and MS
	int megabytes = 0;
	tiff_octaves_init_implicit(v->tg, fgtif, megabytes);
	if (fctif)
		tiff_octaves_init_implicit(v->tc, fctif, megabytes);
	else v->gray_only = 1;//fail("not fctif\n"); //v->tc = NULL;
	v->w = -1;//v->tg->i->w;
	v->h = -1;//v->tg->i->h;
	v->rgbiox = v->rgbioy = 4; // normally untouched

	// load P's RPC
	read_rpc_file_xml(v->r, fgrpc);

	// setup display cache
	v->display = NULL;
	v->fdisplay = NULL;
	v->repaint = 1;
}

static void init_view_no_preview(struct pan_view *v,
		char *fgtif, char *fgrpc, char *fctif, char *fcrpc)
{
	// build preview (use only color, by now)
	//v->preview = load_nice_preview(fgpre, fcpre, &v->pw, &v->ph);
	v->preview = NULL;
	v->pfg = NULL;
	v->pfc = NULL;
	v->gray_only = 0;

	// load P and MS
	int megabytes = 0;
	tiff_octaves_init_implicit(v->tg, fgtif, megabytes);
	if (fctif)
		tiff_octaves_init_implicit(v->tc, fctif, megabytes);
	else v->gray_only = 1;//fail("not fctif\n"); //v->tc = NULL;
	v->w = -1;//v->tg->i->w;
	v->h = -1;//v->tg->i->h;
	v->rgbiox = v->rgbioy = 4; // normally untouched

	// load RPC
	read_rpc_file_xml(v->r, fgrpc);
	if (fcrpc)
		read_rpc_file_xml(v->rc, fcrpc);

	{ // adjust PAN-MSI offsets
		double msi_cx = 2000; //v->tc->i->w / 2.0;
		double msi_cy = 2000; //v->tc->i->h / 2.0;
		double ll[2], ij[2], base_h = 200;
		eval_rpc (ll, v->rc, msi_cx, msi_cy, base_h);
		eval_rpci(ij, v->r , ll[0] , ll[1] , base_h);
		fprintf(stderr, "msi_center = %g %g\n", msi_cx, msi_cy);
		fprintf(stderr, "llh = %g %g %g\n", ll[0], ll[1], base_h);
		fprintf(stderr, "ij  = %g %g\n", ij[0], ij[1]);
		v->rgbiox = 4 * msi_cx - ij[0];
		v->rgbioy = 4 * msi_cy - ij[1];
		fprintf(stderr, "rgbio  = %g %g\n", v->rgbiox, v->rgbioy);
	}

	// setup display cache
	v->display = NULL;
	v->fdisplay = NULL;
	v->repaint = 1;
}

static void init_view_sizes(struct pan_view *v)
{
	void *p = tiff_octaves_getpixel(v->tg, 0, 1, 1);
	if (!p) fail("could not load view %p", (void*)v);
	v->w = v->tg->i->w;
	v->h = v->tg->i->h;
}

static void setup_nominal_pixels_according_to_first_view(struct pan_state *e)
{
	struct pan_view *v = e->view + 0;

	init_view_sizes(v);

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
	double lonfactor = 1.0*cos(latitude);
	double lon_step = -lat_step / lonfactor;

	// fill-in the fields
	e->lat_d = lat_step;
	e->lon_d = lon_step;
	e->lon_0 = center[0];
	e->lat_0 = center[1];
}

static void init_state(struct pan_state *e,
	char *gi[], char *gp[], char *gr[], char *ci[], char *cp[], int n)
{
	// set views
	assert(n < MAX_VIEWS);
	e->nviews = n;
	e->current_view = 0;
	for (int i = 0; i < n; i++)
		init_view(e->view + i,
				gi[i],
				gp?gp[i]:NULL,
				gr[i],
				ci?ci[i]:NULL,
				cp?cp[i]:NULL);

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
	e->force_exact = 0;
	e->image_space = 0;
	e->diff_mode = 0;
	e->show_vertdir = 0;
	e->interpolation_order = 0;
	e->image_rotation_status = 0;
	e->qauto = 0;
	e->log_scale = 0;
	e->show_srtm4 = 0;
	e->so = 0;
	e->srtm4_base = 0;
	msoctaves_instead_of_preview = e->view->tg->noctaves > 4 ||
		(ci && e->view->tc->noctaves > 3);
}

static void init_state_no_preview(struct pan_state *e,
	char *gi[], char *gr[], char *ci[], char *cr[], int n)
{
	// set views
	assert(n < MAX_VIEWS);
	e->nviews = n;
	e->current_view = 0;
	for (int i = 0; i < n; i++)
		init_view_no_preview(e->view + i,
				gi[i],
				gr[i],
				ci?ci[i]:NULL,
				cr?cr[i]:NULL);

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
	e->force_exact = 0;
	e->image_space = 0;
	e->diff_mode = 0;
	e->show_vertdir = 0;
	e->interpolation_order = 0;
	e->image_rotation_status = 0;
	e->qauto = 0;
	e->log_scale = 0;
	e->show_srtm4 = 0;
	e->so = 0;
	e->srtm4_base = 0;
	msoctaves_instead_of_preview = true;e->view->tg->noctaves > 4 ||
		(ci && e->view->tc->noctaves > 3);
}

// state query functions {{{1
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
	if (!msoctaves_instead_of_preview) { // regular case
		if (e->zoom_factor < 0.18) return 0;
		if (e->zoom_factor > 0.5) return 2;
		return 1;
	} else {
		if (e->zoom_factor > 0.5) return 0;
		if (e->zoom_factor > 0.18) return 1;
		if (e->zoom_factor > 0.1) return 2;
		int r = floor( -log2( e->zoom_factor ) );
		fprintf(stderr, "oo z=%g r=%d\n", e->zoom_factor, r);
		return r;
	}
}


// affine approximation of the projection function {{{1
static void window_to_image_exh(double*,struct pan_state*,double,double,double);

static void raster_to_image_exh(double out[2], struct pan_state *e,
		double x, double y, double h);

typedef void raster_to_image_t(double out[2], struct pan_state *e,
		double x, double y, double h);

// compute an affine approximation of a generic raster_to_image function
static void approximate_projection_gen(double P[8],
		struct pan_state *e,
		double x, double y, double h,
		raster_to_image_t *rtoi)
{
	//fprintf(stderr, "gen approximate projection rxyh=%g %g %g\n", x,y,h);
	// finite difference steps
	double eps_R = 0.5; // (in raster steps)
	double eps_H = 1.0; // (in meters)

	// eval the function at a coordinate tetrahedron
	double v0[3], vx[3], vy[3], vh[3];
	rtoi(v0, e, x         , y        , h        );
	rtoi(vx, e, x + eps_R , y        , h        );
	rtoi(vy, e, x         , y + eps_R, h        );
	rtoi(vh, e, x         , y        , h + eps_H);

	// compute forward differences
	double pp = v0[0];
	double qq = v0[1];
	double px = ( vx[0] - v0[0] ) / eps_R;
	double py = ( vy[0] - v0[0] ) / eps_R;
	double ph = ( vh[0] - v0[0] ) / eps_H;
	double qx = ( vx[1] - v0[1] ) / eps_R;
	double qy = ( vy[1] - v0[1] ) / eps_R;
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



// coordinate changes {{{1

// There are four (!) coordinate systems used in this program.
// (I believe that this is the simplest solution, however the code could be
// shortened by identifying the "raster" and the "window" coordinates.)
//
// 1. IMAGE coordinates (p,q), integers in the range [0,40.000] representing the
// position of a pixel in the "P" image.
//
// 2. GEO coordinates (lon,lat), floating point numbers in degrees, e.g., in an
// interval such as [ -137.4782 , -137.4779 ].  These coordinates are the
// geographic position around the site of interest.
//
// 3. RASTER coordinates (x,y), floating point numbers in a small range
// [-20.000, 20.000] representing an isotropic grid over the geographic site,
// approximately at the nominal pixel resolution (may be rotated, scaled, or
// something).
//
// 4. WINDOW coordintes (i,j), integers representing the pixel position in the
// display window

static
void window_to_raster(double xy[2], struct pan_state *e, double i, double j)
{
	xy[0] = e->offset_x + i / e->zoom_factor;
	xy[1] = e->offset_y + j / e->zoom_factor;
}

static
void raster_to_window(double ij[2], struct pan_state *e, double x, double y)
{
	ij[0] = ( x - e->offset_x ) * e->zoom_factor;
	ij[1] = ( y - e->offset_y ) * e->zoom_factor;
}

static
void raster_to_geo(double ll[2], struct pan_state *e, double x, double y)
{
	ll[0] = e->lon_0 + x * e->lon_d;
	ll[1] = e->lat_0 + y * e->lat_d;
}

static
void geo_to_raster(double xy[2], struct pan_state *e, double lon, double lat)
{
	xy[0] = ( lon - e->lon_0 ) / e->lon_d;
	xy[1] = ( lat - e->lat_0 ) / e->lat_d;
}

static
void window_to_image_apm(double out[2], struct pan_state *e, double i, double j)
{
	double xyh[3];
	window_to_raster(xyh, e, i, j);
	double hhh = e->base_h;
	//if (e->srtm4_base) {
	//	double ll[2];
	//	raster_to_geo(ll, e, xyh[0], xyh[1]);
	//	hhh = srtm4(ll[0], ll[1]) + egm96(ll[0], ll[1]);
	//	hhh += e->base_h;
	//}
	xyh[2] = hhh;

	struct pan_view *v = obtain_view(e);
	double p[3];
	apply_projection(p, v->P, xyh);
	out[0] = p[0];
	out[1] = p[1];
}

static
void image_to_window_apm(double ij[2], struct pan_state *e, double p, double q)
{
	struct pan_view *v = obtain_view(e);
	double *P = v->P;
	double a = P[0];
	double b = P[1];
	double c = P[4];
	double d = P[5];
	double det = a * d - b * c;
	double x = p - P[2]*e->base_h - P[3];
	double y = q - P[6]*e->base_h - P[7];
	double raster[2];
	raster[0] = ( d * x - b * y ) / det;
	raster[1] = ( a * y - c * y ) / det;
	raster_to_window(ij, e, raster[0], raster[1]);
}


// from the geographic raster to the image domain, using the RPC functions
static void raster_to_image_exh(double out[2], struct pan_state *e,
		double x, double y, double h)
{
	double lonlat[2];
	raster_to_geo(lonlat, e, x, y);

	struct pan_view *v = obtain_view(e);
	double p[3];
	eval_rpci(p, v->r, lonlat[0], lonlat[1], h);
	out[0] = p[0];
	out[1] = p[1];
}

static void vertical_direction(double vertical[2],
		struct pan_view *v, double p, double q, double h)
{
	p = v->w / 2;
	q = v->h / 2;
	double lonlat[3];
	eval_rpc(lonlat, v->r, p, q, h);


	double eps_H = 100.0; // in meters
	double p_bot[3], p_top[3];
	eval_rpci(p_bot, v->r, lonlat[0], lonlat[1], h        );
	eval_rpci(p_top, v->r, lonlat[0], lonlat[1], h + eps_H);
	vertical[0] = ( p_top[0] - p_bot[0] ) / eps_H;
	vertical[1] = ( p_top[1] - p_bot[1] ) / eps_H;
}

// from the image raster to the image coordinates
// (typically, a translation along the vertical direction, proportional to h)
//
// This function is slow, intended to be approximated later by an affine map
static void raster_to_image_raw(double out[2], struct pan_state *e,
		double rx, double ry, double h)
{
	struct pan_view *v = obtain_view(e);

	// image point corresponding to raster x,y
	// (HERE identity map)
	//double x = -rx;
	//double y = -ry;

	double vertical[2];
	vertical_direction(vertical, v, rx, ry, h);

	double alpha = atan2(vertical[1], vertical[0]) + M_PI/2;
	if (e->image_rotation_status == 0) alpha = 0;
	if (e->image_rotation_status == 1) alpha = -M_PI;
	double x = cos(alpha) * rx - sin(alpha) * ry;
	double y = sin(alpha) * rx + cos(alpha) * ry;


	out[0] = x + h*vertical[0];
	out[1] = y + h*vertical[1];
}

static void image_to_raster_raw(double out[2], struct pan_state *e,
		double p, double q, double h)
{
	struct pan_view *v = obtain_view(e);

	double vertical[2];
	vertical_direction(vertical, v, p, q, h);

	double x = p - h * vertical[0];
	double y = q - h * vertical[1];

	// (HERE inverse identity map)
	//double alpha = 0.1;
	double alpha = atan2(vertical[1], vertical[0]) + M_PI/2;
	if (e->image_rotation_status == 0) alpha = 0;
	if (e->image_rotation_status == 1) alpha = -M_PI;
	out[0] =  cos(alpha) * x + sin(alpha) * y;
	out[1] = -sin(alpha) * x + cos(alpha) * y;
}

static
void window_to_image_exh(double out[2], struct pan_state *e,
		double i, double j, double h)
{
	double xy[2];
	window_to_raster(xy, e, i, j);
	raster_to_image_exh(out, e, xy[0], xy[1], h);
}

static
void window_to_image_ex(double out[2], struct pan_state *e, double i, double j)
{
	window_to_image_exh(out, e, i, j, e->base_h);
}

static
void image_to_window_ex(double ij[2], struct pan_state *e, double p, double q)
{
	struct pan_view *v = obtain_view(e);
	double lonlat[3], xy[2];
	eval_rpc(lonlat, v->r, p, q, e->base_h);
	geo_to_raster(xy, e, lonlat[0], lonlat[1]);
	raster_to_window(ij, e, xy[0], xy[1]);
}


// pixel interpolation {{{1

// constant, bilinear, bicubic {{{2
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

static void preview_at_bil(float *out, struct pan_view *v, double p, double q)
{
	p -= 0.5;
	q -= 0.5;
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
		out[l] = r;
	}
}

static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}


static void preview_at_bic(float *out, struct pan_view *v, double x, double y)
{
	x -= 1.5;
	y -= 1.5;
	int ix = floor(x);
	int iy = floor(y);
	int w = v->pw;
	int h = v->ph;
	unsigned char *img = v->preview;
	for (int l = 0; l < 3; l++) {
		float c[4][4];
		for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = rgb_getsamplec(img, w, h, ix + i, iy + j, l);
		out[l] = bicubic_interpolation_cell(c, x - ix, y - iy);
	}
}

// RGBI and MS combination stuff {{{2
//static float getgreenf(float c[4])
//{
//	return c[1] * 0.6 + c[3] * 0.2;
//}

static void rgbi_to_rgb_inplace(float c[4])
{
	float t = c[0];
	c[0] = c[2];
	c[2] = t;
	c[1] = 0.8*c[1] + 0.2*c[3];
	//c[0] = 1.00 * c[0]  +  0.05 * c[3];
	//c[1] = 0.60 * c[1]  +  0.20 * c[3];
	//c[2] = 1.30 * c[2]  -  0.20 * c[3];
}

static void rgbi_to_rgb_inplace_n(float *c, int n)
{
	float t;
	switch (n)
	{
	case 4:
		t = c[0];
		c[0] = c[2];
		c[2] = t;
		c[1] = 0.8*c[1] + 0.2*c[3];
		break;
	case 8:
		c[0] = c[4];
		t = c[1];
		c[1] = c[2];
		c[2] = t;
		break;
	default: fail("caca");
	}
	//c[0] = 1.00 * c[0]  +  0.05 * c[3];
	//c[1] = 0.60 * c[1]  +  0.20 * c[3];
	//c[2] = 1.30 * c[2]  -  0.20 * c[3];
}

static
void pixel_from_preview(float *r, struct pan_view *v, double p, double q, int i)
{
	if (!v->preview)
		v->preview = load_nice_preview(v->pfg, v->pfc, &v->pw, &v->ph);
	int ip = p * v->pw / v->tg->i->w;
	int iq = q * v->ph / v->tg->i->h;
	float fp = p * v->pw / v->tg->i->w;
	float fq = q * v->ph / v->tg->i->h;
	if (insideP(v->pw, v->ph, ip, iq)) {
		if (i == 2)      preview_at_bil(r, v, fp, fq);
		else if (i == 3) preview_at_bic(r, v, fp, fq);
		else for (int l = 0; l < 3; l++)
			r[l] = rgb_getsamplec(v->preview, v->pw, v->ph,
					ip, iq, l);
	} else {
		// bright green
		r[0] = 0;
		r[1] = 255;
		r[2] = 0;
	}
}

static int tiffo_getpixel_float_raw(float *r, struct tiff_octaves *t,
		int o, int i, int j)
{
	void *p = tiff_octaves_getpixel(t, o, i, j);
	convert_pixel_to_float(r, t->i, p);
	return t->i->spp;
}

static int tiffo_getpixel_float_1(float *r, struct tiff_octaves *t,
		int o, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= t->i[o].w) i = t->i[o].w;
	if (j >= t->i[o].h) i = t->i[o].h;
	return tiffo_getpixel_float_raw(r, t, o, i, j);
}

static int tiffo_getpixel_float_bilinear(float *r, struct tiff_octaves *t,
		int o, float p, float q)
{
	p -= 0.5;
	q -= 0.5;
	int ip = p;
	int iq = q;
	int pd = t->i->spp;
	float v[4][pd];
	tiffo_getpixel_float_1(v[0], t, o, ip  , iq  );
	tiffo_getpixel_float_1(v[1], t, o, ip+1, iq  );
	tiffo_getpixel_float_1(v[2], t, o, ip  , iq+1);
	tiffo_getpixel_float_1(v[3], t, o, ip+1, iq+1);
	for (int l = 0; l < pd; l++)
		r[l] = evaluate_bilinear_cell(v[0][l], v[1][l],
						v[2][l], v[3][l], p-ip, q-iq);
	return pd;
}

static int tiffo_getpixel_float_bicubic(float *r, struct tiff_octaves *t,
		int o, float x, float y)
{
	x -= 1.5;
	y -= 1.5;
	int ix = floor(x);
	int iy = floor(y);
	int pd = t->i->spp;
	float c[4][4][pd];
	for (int j = 0; j < 4; j++)
	for (int i = 0; i < 4; i++)
		tiffo_getpixel_float_1(c[i][j], t, o, ix + i, iy + j);
	for (int l = 0; l < pd; l++) {
		float cl[4][4];
		for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			cl[i][j] = c[i][j][l];
		r[l] = bicubic_interpolation_cell(cl, x - ix, y - iy);
	}
	return pd;
}

static void pixel_from_ms(float *out, struct pan_view *v, double p, double q,
		int i)
{
	float c[20];
	float ipc = (p + v->rgbiox) / 4;
	float iqc = (q + v->rgbioy) / 4;
	if (i == 2)      tiffo_getpixel_float_bilinear(c, v->tc, 0, ipc, iqc);
	else if (i == 3) tiffo_getpixel_float_bicubic(c, v->tc, 0, ipc, iqc);
	else             tiffo_getpixel_float_1(c, v->tc, 0, ipc, iqc);
	rgbi_to_rgb_inplace_n(c, v->tc->i->spp);
	float g = (c[0] + c[1] + c[2] + c[3]) / 4;
	float nc = 4*hypot(c[0], hypot(c[1], c[2]));
	out[0] = c[0] * g / nc;
	out[1] = c[1] * g / nc;
	out[2] = c[2] * g / nc;
}

static void pixel_from_mso(float *out, struct pan_view *v, double p, double q,
		int i, int o)
{
	////
	//out[0] = 0;
	//out[1] = 50;
	//out[2] = 200;
	int ofac = 1 << (o);
	float c[20];
	float ipc = (p + v->rgbiox) / ofac;
	float iqc = (q + v->rgbioy) / ofac;
	if (i == 2)      tiffo_getpixel_float_bilinear(c, v->tc, o-2, ipc, iqc);
	else if (i == 3) tiffo_getpixel_float_bicubic(c, v->tc, o-2, ipc, iqc);
	else             tiffo_getpixel_float_1(c, v->tc, o-2, ipc, iqc);
	//rgbi_to_rgb_inplace(c);
	rgbi_to_rgb_inplace_n(c, v->tc->i->spp);
	float g = (c[0] + c[1] + c[2] + c[3]) / 4;
	float nc = 4*hypot(c[0], hypot(c[1], c[2]));
	out[0] = c[0] * g / nc;
	out[1] = c[1] * g / nc;
	out[2] = c[2] * g / nc;
}

static void pixel_from_go(float *out, struct pan_view *v, double p, double q,
		int i, int o)
{
	int ofac = 1 << (o);
	float c[20];
	float ipc = p / ofac;
	float iqc = q / ofac;
	if (i == 2)      tiffo_getpixel_float_bilinear(c, v->tg, o, ipc, iqc);
	else if (i == 3) tiffo_getpixel_float_bicubic(c, v->tg, o, ipc, iqc);
	else             tiffo_getpixel_float_1(c, v->tg, o, ipc, iqc);
	out[0] = out[1] = out[2] = c[0];
}

static void pixel_from_pms(float *out, struct pan_view *v, double p, double q,
		int i)
{
	float pc = (p + v->rgbiox) / 4;
	float qc = (q + v->rgbioy) / 4;
	float g, c[20];
	if (i == 2)      tiffo_getpixel_float_bilinear(&g, v->tg, 0, p , q );
	else if (i == 3) tiffo_getpixel_float_bicubic (&g, v->tg, 0, p , q );
	else             tiffo_getpixel_float_1       (&g, v->tg, 0, p , q );
	if (i == 2)      tiffo_getpixel_float_bilinear(c , v->tc, 0, pc, qc);
	else if (i == 3) tiffo_getpixel_float_bicubic (c , v->tc, 0, pc, qc);
	else             tiffo_getpixel_float_1       (c , v->tc, 0, pc, qc);
	//rgbi_to_rgb_inplace(c);
	rgbi_to_rgb_inplace_n(c, v->tc->i->spp);
	float nc = 4*hypot(c[0], hypot(c[1], c[2]));
	out[0] = c[0] * g / nc;
	out[1] = c[1] * g / nc;
	out[2] = c[2] * g / nc;
}

// generic "pixel" (includes all the previous cases) {{{2

// evaluate the value a position (p,q) in image coordinates
static
void pixel(float *out, struct pan_view *v, double p, double q, int o, int i)
{
	// TODO: remove the conditionals inside this function by putting
	// all the version inside a table of pointers to functions
	if (p < 0 || q < 0 || p >= v->w || q >= v->h) {
		// pink
		out[0] = 200;
		out[1] = 50;
		out[2] = 100;
		return;
	}
	if (v->gray_only)
		pixel_from_go(out, v, p, q, i, o);
	else if (msoctaves_instead_of_preview) { // ms-octaves case
		if (o == 0 || o == 1)
			pixel_from_pms (out, v, p, q, i);
		else if (o == 2)
			pixel_from_ms  (out, v, p, q, i);
		else
			pixel_from_mso (out, v, p, q, i, o);
	} else { // regular case
		if (o == 0)      pixel_from_preview (out, v, p, q, i);
		else if (o == 1) pixel_from_ms      (out, v, p, q, i);
		else if (o == 2) pixel_from_pms     (out, v, p, q, i);
		else exit(fprintf(stderr,"ERROR: bad octave %d\n", o));
	}
}

// contrast changes {{{1


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

static void inplace_rgb_span(float *x, int w, int h, double a)
{
	int n = w*h;
	float min[3] = {INFINITY, INFINITY, INFINITY}, max[3] = {0, 0, 0};
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
	{
		min[l] = fmin(min[l], x[3*i+l]);
		max[l] = fmax(max[l], x[3*i+l]);
	}
	fprintf(stderr, "qauto %g %g %g | %g %g %g\n",
			min[0], min[1], min[2], max[0], max[1], max[2]);
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
	{
		float xx = ( x[3*i+l] - min[l] ) / ( max[l] - min[l] );
		x[3*i+l] = 255.0 * pow( xx, a );
	}
}

static void inplace_rgb_span2(float *x, int w, int h, double a)
{
	double avg[3] = {0, 0, 0}, var[3] = {0, 0, 0};
	int mw = w/4, mh = h/4;
	int n = (w-2*mw)*(h-2*mh);
	for (int j = mh; j < h - mh ; j++)
	for (int i = mw; i < w - mw ; i++)
	for (int l = 0; l < 3; l++)
		avg[l] += x[3*(j*w+i)+l]/n;
	for (int j = mh; j < h - mh ; j++)
	for (int i = mw; i < w - mw ; i++)
	for (int l = 0; l < 3; l++)
	{
		double q = x[3*(j*w+i)+l] - avg[l];
		var[l] += q * q / n;
	}
	fprintf(stderr, "span2 avg = %g %g %g | dev = %g %g %g\n",
				avg[0], avg[1], avg[2],
				sqrt(var[0]), sqrt(var[1]), sqrt(var[2]));
	if (isfinite(*avg) && isfinite(*var))
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < 3; l++)
		x[3*(j*w+i)+l] = 127 + a*(x[3*(j*w+i)+l] - avg[l])/sqrt(var[l]);
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void inplace_rgb_span3(float *x, int w, int h, double a)
{
	int n = w*h;
	int nnn[3] = {0, 0, 0};
	float *tmp[3];
	for (int l = 0; l < 3; l++) tmp[l] = xmalloc(n*sizeof(float));
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
		if (isfinite(x[3*i+l]))
			tmp[l][nnn[l]++] = x[3*i+l];
	for (int l = 0; l < 3; l++)
		qsort(tmp[l], n, sizeof(float), compare_floats);
	float med[3], iqd[3];
	for (int l = 0; l < 3; l++) med[l] = tmp[l][n/2];
	for (int l = 0; l < 3; l++) iqd[l] = tmp[l][3*n/4] - tmp[l][1*n/4];
	fprintf(stderr, "qauto med %g %g %g | iqd %g %g %g\n",
			med[0], med[1], med[2], iqd[0], iqd[1], iqd[2]);
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
	{
		float xx = ( x[3*i+l] - med[l] ) / iqd[l];
		x[3*i+l] = 127 + a * (x[3*i+l] - med[l] ) / iqd[l];
	}
	for (int l = 0; l < 3; l++) free(tmp[l]);
}


// pan_repaint {{{1
//
static void update_local_projection(struct pan_state *e,
		double ix, double iy, double h)
{
	struct pan_view *v = obtain_view(e);

	double xy[2];
	window_to_raster(xy, e, ix, iy);
	double x = xy[0];
	double y = xy[1];

	approximate_projection_gen(v->P, e, x, y, h, e->image_space ?
			raster_to_image_raw : raster_to_image_exh);
}


// dump the image acording to the state of the viewport
static void pan_repaint(struct pan_state *e, int w, int h)
{
	struct pan_view *v = obtain_view(e);

	// setup the display buffer
	if (!v->display || v->dw != w || v->dh != h) {
		if (v->display)
			free(v->display);
		v->display = xmalloc(w * h * 3);
		v->fdisplay = xmalloc(w * h * 3 * sizeof*v->fdisplay);
		v->dw = w;
		v->dh = h;
		v->repaint = 1;
	}
	if (!v->repaint) return; // if no repaint requested, return
	v->repaint = 0;

	double dh = 0;
	if (e->srtm4_base) {
		abort();
		double xy[2], ll[2];
		window_to_raster(xy, e, w/2, h/2);
		raster_to_geo(ll, e, xy[0], xy[1]);
		dh = srtm4o(ll[0], ll[1], 0) + egm96(ll[0], ll[1]);
		e->base_h = dh;
	}
	update_local_projection(e, w/2, h/2, e->base_h);

	static void (*win_to_img)(double p[2],struct pan_state*,double,double);
	win_to_img = window_to_image_apm;
	if (e->image_space)      win_to_img = window_to_image_apm;
	else if (e->force_exact) win_to_img = window_to_image_ex;

	int o = obtain_octave(e);
	int interp = e->interpolation_order;

//#pragma omp parallel for
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!e->image_space && e->show_srtm4) {
			double xy[3], ll[3];
			window_to_raster(xy, e, i, j);
			raster_to_geo(ll, e, xy[0], xy[1]);
			int so = -8 + lrint(-log2(e->zoom_factor));
			if (so < 0) so = 0;
			if (!i && !j)
			{
				fprintf(stderr, "z=%g, so=%d ofac=%d\n",
						e->zoom_factor, so, 1 << so);
				so = -so;
			}
			double hhh = srtm4o(ll[0], ll[1], so) + egm96(ll[0], ll[1]);
			float *cc = v->fdisplay + 3 * (j * v->dw + i);
			for (int l = 0; l < 3; l++)
				cc[l] = e->a * hhh + e->b;
		} else if (!e->diff_mode) {
			double p[2];
			win_to_img(p, e, i, j);
			float c[3];
			pixel(c, v, p[0], p[1], o, interp);
			float *cc = v->fdisplay + 3 * (j * v->dw + i);
			if (!e->log_scale)
				for (int l = 0; l < 3; l++)
					cc[l] = e->a * c[l] + e->b;
			else
				for (int l = 0; l < 3; l++)
					cc[l] = e->a*(32*log2(c[l]+1)-1) + e->b;
		} else { // diff mode
			int va = e->current_view;
			int vb = good_modulus(e->current_view+1,e->nviews);
			double p[2], q[2];
			e->current_view = va; win_to_img(p, e, i, j);
			e->current_view = vb; win_to_img(q, e, i, j);
			float ca[3], cb[3];
			pixel(ca, e->view + va, p[0], p[1], o, interp);
			pixel(cb, e->view + vb, q[0], q[1], o, interp);
			float *cc = v->fdisplay + 3 * (j * v->dw + i);
			for (int l = 0; l < 3; l++) {
				float g = 2*e->a * (ca[l] - cb[l]) + 127;
				cc[l] = g;
			}
			e->current_view = va;
		}
	}
	if (e->qauto == 1) inplace_rgb_span(v->fdisplay, v->dw, v->dh, 1/e->a);
	if (e->qauto == 2)inplace_rgb_span2(v->fdisplay, v->dw, v->dh, 40*e->a);
	if (e->qauto == 3)inplace_rgb_span3(v->fdisplay, v->dw, v->dh, 40*e->a);
	for (int i = 0; i < v->dw * v->dh * 3; i++)
		v->display[i] = float_to_uint8(v->fdisplay[i]);
}

// CALLBACK: pan_exposer {{{1

static void overlay_vertdir(struct FTR *f, int i, int j)
{
	struct pan_state *e = f->userdata;
	fprintf(stderr, "vertdir %d %d\n", i, j);

	for (double height = 0; height < 100; height += 1)
	{
		double pq[2], xy[2];
		window_to_image_exh(pq, e, i, j, e->base_h + height);
		image_to_window_ex(xy, e, pq[0], pq[1]);

		if (insideP(f->w, f->h, xy[0], xy[1])) {
			f->rgb[3*(f->w*lrint(xy[1]) + lrint(xy[0])) + 0] = 255;
			f->rgb[3*(f->w*lrint(xy[1]) + lrint(xy[0])) + 1] = 0;
			f->rgb[3*(f->w*lrint(xy[1]) + lrint(xy[0])) + 2] = 0;
		}
	}
}

static void pan_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	(void)unused_x; (void)unused_y;
	//fprintf(stderr, "\n\nexpose %d %d\n", b, m);
	struct pan_state *e = f->userdata;
	struct pan_view *v = obtain_view(e);

	pan_repaint(e, f->w, f->h);

	// copy the requested view into the display
	assert(f->w == v->dw);
	assert(f->h == v->dh);
	memcpy(f->rgb, v->display, v->dw * v->dh * 3);

	// overlay any requested lines
	if (e->show_vertdir)
		overlay_vertdir(f, e->vdx, e->vdy);
}

static void request_repaints(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	for (int i = 0; i < e->nviews; i++)
		e->view[i].repaint = 1;
	f->changed = 1;
}

// Actions {{{1

// TODO: edit all actions so that they act directly on the "pan_state", not on
// the "FTR".  This is necessary so that actions can be directly called even
// when there is no interface (e.g., for scripted usage).

static void action_offset_viewport(struct FTR *f, double dx, double dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->zoom_factor;
	e->offset_y -= dy/e->zoom_factor;

	request_repaints(f);
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

	request_repaints(f);
}

static void action_multiply_contrast(struct FTR *f, double fac)
{
	struct pan_state *e = f->userdata;
	e->a *= fac;
	fprintf(stderr, "\t constrast factor changed to %g\n", e->a);

	request_repaints(f);
}

static void action_offset_base_h(struct FTR *f, double d)
{
	struct pan_state *e = f->userdata;
	e->base_h += d;
	fprintf(stderr, "base_h = %g\n", e->base_h);

	request_repaints(f);
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

	request_repaints(f);
}

static void action_toggle_srtm4(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->show_srtm4 = !e->show_srtm4;
	request_repaints(f);
}

static void action_shift_so(struct FTR *f, int dso)
{
	struct pan_state *e = f->userdata;
	e->so += dso;
	request_repaints(f);
}

static void action_toggle_srtm4_base(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->srtm4_base = !e->srtm4_base;
	if (e->srtm4_base)
		fprintf(stderr, "use SRTM4 base_h\n");
	else
		fprintf(stderr, "use base_h = %g\n", e->base_h);
	request_repaints(f);
}

static void action_toggle_diff_mode(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->diff_mode = !e->diff_mode;
	request_repaints(f);
}

static void action_toggle_log_scale(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->log_scale = !e->log_scale;
	fprintf(stderr, "log scale = %s\n", e->log_scale?"YES":"NO");
	request_repaints(f);
}

static void action_toggle_vertdir(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;
	e->show_vertdir = !e->show_vertdir;
	e->vdx = x;
	e->vdy = y;
	f->changed = 1;
}

static void action_toggle_image_space(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	fprintf(stderr, "doing hacky things\n");

	// assure that the center of the window is not moved
	e->image_space = !e->image_space;
	if (e->image_space) {
		update_local_projection(e, f->w/2, f->h/2, e->base_h);
		double c[2] = {x, y}, p[2], q[2];
		window_to_image_ex(p, e, c[0], c[1]);
		image_to_window_apm(q, e, p[0], p[1]);
		action_offset_viewport(f, c[0] - q[0], c[1] - q[1]);
		fprintf(stderr, "going window to image around %d %d\n", x, y);
		fprintf(stderr, "\tp = %g %g\n", p[0], p[1]);
		fprintf(stderr, "\tq = %g %g\n", q[0], q[1]);
		fprintf(stderr, "\toff = %g %g\n", c[0]-q[0], c[1]-q[1]);
	} else {
		double c[2] = {x, y}, p[2], q[2];
		window_to_image_apm(p, e, c[0], c[1]);
		image_to_window_ex(q, e, p[0], p[1]);
		action_offset_viewport(f, c[0] - q[0], c[1] - q[1]);
		fprintf(stderr, "going image to window around %d %d\n", x, y);
		fprintf(stderr, "\tp = %g %g\n", p[0], p[1]);
		fprintf(stderr, "\tq = %g %g\n", q[0], q[1]);
		fprintf(stderr, "\toff = %g %g\n", c[0]-q[0], c[1]-q[1]);
	}
	if (e->image_space)
		fprintf(stderr, "base grid = IMAGE\n");
	else
		fprintf(stderr, "base grid = GEOGRAPHIC\n");
	request_repaints(f);
}

static void action_cycle_rotation_status(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->image_rotation_status = ( e->image_rotation_status + 1 ) % 3;
	fprintf(stderr, "rotation status = ");
	if (e->image_rotation_status == 0) fprintf(stderr, "UP\n");
	if (e->image_rotation_status == 1) fprintf(stderr, "DOWN\n");
	if (e->image_rotation_status == 2) fprintf(stderr, "VERTICAL\n");
	// Missing: NORTH, SOUTH, SUN_UP, SUN_DOWN
	// More seriously: this should be independent of the "raster type"
	request_repaints(f);
}

static void action_cycle_contrast(struct FTR *f, int dir)
{
	struct pan_state *e = f->userdata;
	e->qauto = ( e->qauto + dir ) % 4;
	fprintf(stderr, "contrast change = ");
	if (e->qauto == 0) fprintf(stderr, "NONE\n");
	if (e->qauto == 1) fprintf(stderr, "MIN-MAX\n");
	if (e->qauto == 2) fprintf(stderr, "AVG-STD\n");
	if (e->qauto == 3) fprintf(stderr, "MED-IQD\n");
	request_repaints(f);
}

// This function is only called when "ignoring" the RPC functions.  Even when
// we ignore the RPCs, we still want to be able to flip the images while
// keeping the center of the image at the correct offset; to compute this
// offset we need the RPC functions.
static void reposition_in_image_space(struct FTR *f, int a, int b, int x, int y)
{
	if (x < 0 || y < 0 || x >= f->w || y >= f->h) { x = f->w/2; y = f->h/2;}
	struct pan_state *e = f->userdata; assert(e->image_space);
	int save_view = e->current_view;
	struct pan_view *va = e->view + a;
	struct pan_view *vb = e->view + b;
	double rc[2], imgac[2], imgbc[2], rdc[2], dc[2], lonlat[2];
	window_to_raster(rc, e, x, y);
	e->current_view = a;
	raster_to_image_raw(imgac, e, rc[0], rc[1], e->base_h);
	eval_rpc(lonlat, va->r, imgac[0], imgac[1], e->base_h);
	eval_rpci(imgbc, vb->r, lonlat[0], lonlat[1], e->base_h);
	e->current_view = b;
	image_to_raster_raw(rdc, e, imgbc[0], imgbc[1], e->base_h);
	raster_to_window(dc, e, rdc[0], rdc[1]);
	double offx = x - dc[0];
	double offy = y - dc[1];
	e->offset_x -= offx/e->zoom_factor;
	e->offset_y -= offy/e->zoom_factor;
	e->current_view = save_view;
}

static void action_select_view(struct FTR *f, int i, int x, int y)
{
	struct pan_state *e = f->userdata;
	if (i >= 0 && i < e->nviews)
	{
		fprintf(stderr, "selecting view %d\n", i);
		init_view_sizes(e->view + i);
		if (e->image_space)
			reposition_in_image_space(f, e->current_view, i, x, y);
		e->current_view = i;
		f->changed = 1;
	}
}

static void action_select_interpolator(struct FTR *f, int k)
{
	struct pan_state *e = f->userdata;
	if (k >= 0 || k < 3 || k == 4)
	       	e->interpolation_order = k;
	request_repaints(f);
}

static void action_offset_rgbi(struct FTR *f, double dx, double dy)
{
	struct pan_state *e = f->userdata;
	struct pan_view *v = obtain_view(e);
	v->rgbiox += dx;
	v->rgbioy += dy;
	fprintf(stderr, "rgbi offset chanted to %g %g\n", v->rgbiox, v->rgbioy);
	request_repaints(f);
}

static void action_cycle_view(struct FTR *f, int d, int x, int y)
{
	struct pan_state *e = f->userdata;
	int new = good_modulus(e->current_view + d, e->nviews);
	action_select_view(f, new, x, y);
}

static void action_dump_view(struct FTR *f)
{
	static int dump_view_counter = 0;
	int pid = getpid();
	char fnamei[FILENAME_MAX], fnamef[FILENAME_MAX];
	char *fmti = "/tmp/rpcflip_dump_%d_%d.png";
	char *fmtf = "/tmp/rpcflip_dump_%d_%d.tiff";
	snprintf(fnamei, FILENAME_MAX, fmti, pid, dump_view_counter);
	snprintf(fnamef, FILENAME_MAX, fmtf, pid, dump_view_counter);
	struct pan_state *e = f->userdata;
	struct pan_view  *v = obtain_view(e);
	iio_write_image_uint8_vec(fnamei, v->display, v->dw, v->dh, 3);
	iio_write_image_float_vec(fnamef, v->fdisplay, v->dw, v->dh, 3);
	fprintf(stderr, "dumped rgb_8 view to file \"%s\"\n", fnamei);
	fprintf(stderr, "dumped rgb_f view to file \"%s\"\n", fnamef);
	dump_view_counter += 1;
}

static void action_dump_collection(struct FTR *f)
{
	static int collection_counter = 0;
	struct pan_state *e = f->userdata;
	for (int i = 0; i < e->nviews; i++)
	{
		action_select_view(f, i, f->w/2, f->h/2);
		pan_repaint(e, f->w, f->h);
		int pid = getpid();
		char fname[FILENAME_MAX];
		char *fmt = "/tmp/rpcflip_collection_%d_%d_%d.png";
		snprintf(fname, FILENAME_MAX, fmt, pid, collection_counter, i);
		struct pan_view  *v = obtain_view(e);
		fprintf(stderr, "going to dum rgb_8 view(%d %d) to file \"%s\"\n", v->dw, v->dh, fname);
		iio_write_image_uint8_vec(fname, v->display, v->dw, v->dh, 3);
	}
	collection_counter += 1;
}

static void action_dump_raw(struct FTR *f)
{
	static int raw_view_counter = 0;
	int pid = getpid();
	char fname_pan[FILENAME_MAX], fname_msi[FILENAME_MAX];
	char *fmt_pan = "/tmp/rpcflip_raw_pan_%d_%d.tiff";
	char *fmt_msi = "/tmp/rpcflip_raw_msi_%d_%d.tiff";
	snprintf(fname_pan, FILENAME_MAX, fmt_pan, pid, raw_view_counter);
	snprintf(fname_msi, FILENAME_MAX, fmt_msi, pid, raw_view_counter);
	struct pan_state *e = f->userdata;
	struct pan_view  *v = obtain_view(e);
	float *buf_pan = xmalloc(1 * v->dw * v->dh * sizeof(float));
	float *buf_msi = xmalloc(4 * v->dw * v->dh * sizeof(float));
	for (int j = 0; j < v->dh; j++)
	for (int i = 0; i < v->dw; i++)
	{
		// basically, redo the function "pan_repaint"
		static
		void (*win_to_img)(double p[2],struct pan_state*,double,double);
		win_to_img = window_to_image_apm;
		if (e->image_space)      win_to_img = window_to_image_apm;
		else if (e->force_exact) win_to_img = window_to_image_ex;
		int o = obtain_octave(e);
		int interp = e->interpolation_order;
		if (!e->image_space) continue;
		double p[2];
		win_to_img(p, e, i, j);
		if (p[0]<0 || p[1]<0 || p[0]>=v->w || p[1]>=v->h)
			continue;
		float pc = (p[0] + v->rgbiox) / 4;
		float qc = (p[1] + v->rgbioy) / 4;
		float g = 0, c[4] = {0, 0, 0, 0};
		if (i == 2) tiffo_getpixel_float_bilinear(&g,v->tg,0,p[0],p[1]);
		else if(i==3)tiffo_getpixel_float_bicubic(&g,v->tg,0,p[0],p[1]);
		else tiffo_getpixel_float_1              (&g,v->tg,0,p[0],p[1]);
		if (!v->gray_only) {
		if (i == 2) tiffo_getpixel_float_bilinear(c,v->tc,0,pc,qc);
		else if(i==3)tiffo_getpixel_float_bicubic(c,v->tc,0,pc,qc);
		else       tiffo_getpixel_float_1        (c,v->tc,0,pc,qc);
		}
		buf_pan[i+j*v->dw] = g;
		for (int l = 0; l < 4; l++)
			buf_msi[(i+j*v->dw)*4+l] = c[l];
	}
	iio_write_image_float_vec(fname_pan, buf_pan, v->dw, v->dh, 1);
	if (!v->gray_only)
	iio_write_image_float_vec(fname_msi, buf_msi, v->dw, v->dh, 4);
	fprintf(stderr, "dumped raw pan to file \"%s\"\n", fname_pan);
	if (!v->gray_only)
	fprintf(stderr, "dumped raw msi to file \"%s\"\n", fname_msi);
	raw_view_counter += 1;
}

static void action_dump_raw_fancy(struct FTR *f)
{
	static int raw_view_counter = 0;
	int pid = getpid();
	char fname_pan[FILENAME_MAX], fname_msi[FILENAME_MAX];
	char *fmt_pan = "/tmp/rpcflip_fraw_pan_%d_%d.tiff";
	char *fmt_msi = "/tmp/rpcflip_fraw_msi_%d_%d.tiff";
	snprintf(fname_pan, FILENAME_MAX, fmt_pan, pid, raw_view_counter);
	snprintf(fname_msi, FILENAME_MAX, fmt_msi, pid, raw_view_counter);
	struct pan_state *e = f->userdata;
	struct pan_view  *v = obtain_view(e);
	int pd = v->tc->i->spp;
	float *buf_pan = xmalloc( 1 * v->dw * v->dh * sizeof(float));
	float *buf_msi = xmalloc(pd * v->dw * v->dh * sizeof(float));
	for (int j = 0; j < v->dh; j++)
	for (int i = 0; i < v->dw; i++)
	{
		// basically, redo the function "pan_repaint"
		static
		void (*win_to_img)(double p[2],struct pan_state*,double,double);
		win_to_img = window_to_image_apm;
		if (e->image_space)      win_to_img = window_to_image_apm;
		else if (e->force_exact) win_to_img = window_to_image_ex;
		int o = obtain_octave(e);
		int io = e->interpolation_order;
		//if (!e->image_space) continue;
		double p[2];
		win_to_img(p, e, i, j);
		if (p[0]<0 || p[1]<0 || p[0]>=v->w || p[1]>=v->h)
			continue;
		float pc = (p[0] + v->rgbiox) / 4;
		float qc = (p[1] + v->rgbioy) / 4;
		float g = i+0.1*j, c[pd];
		if (io == 2)tiffo_getpixel_float_bilinear(&g,v->tg,0,*p,p[1]);
		else if(io==3)tiffo_getpixel_float_bicubic(&g,v->tg,0,*p,p[1]);
		else tiffo_getpixel_float_1              (&g,v->tg,0,p[0],p[1]);
		if (!v->gray_only) {
		if (io == 2) tiffo_getpixel_float_bilinear(c,v->tc,0,pc,qc);
		else if(io==3)tiffo_getpixel_float_bicubic(c,v->tc,0,pc,qc);
		else       tiffo_getpixel_float_1        (c,v->tc,0,pc,qc);
		}
		buf_pan[i+j*v->dw] = g;
		for (int l = 0; l < pd; l++)
			buf_msi[(i+j*v->dw)*pd+l] = c[l];
	}
	iio_write_image_float_vec(fname_pan, buf_pan, v->dw, v->dh,  1);
	if (!v->gray_only)
	iio_write_image_float_vec(fname_msi, buf_msi, v->dw, v->dh, pd);
	fprintf(stderr, "dumped raw pan to file \"%s\"\n", fname_pan);
	if (!v->gray_only)
	fprintf(stderr, "dumped raw msi (pd=%d) to file \"%s\"\n",pd,fname_msi);
	raw_view_counter += 1;
}

static void action_dump_raw_collection_fancy(struct FTR *f)
{
	static int raw_collection_counter = 0;
	struct pan_state *e = f->userdata;
	int pid = getpid();
	int F=FILENAME_MAX;
	char fname_pan[F], fname_msi[F];
	char *fmt_pan = "/tmp/rpcflip_cfraw_pan_%d_%d_%d.tiff";
	char *fmt_msi = "/tmp/rpcflip_cfraw_msi_%d_%d_%d.tiff";

	for (int i = 0; i < e->nviews; i++)
	{
		action_select_view(f, i, f->w/2, f->h/2);
		pan_repaint(e, f->w, f->h);
		snprintf(fname_pan, F, fmt_pan, pid, raw_collection_counter, i);
		snprintf(fname_msi, F, fmt_msi, pid, raw_collection_counter, i);
		struct pan_view  *v = obtain_view(e);

		int pd = v->tc->i->spp;
		float *buf_pan = xmalloc( 1 * v->dw * v->dh * sizeof(float));
		float *buf_msi = xmalloc(pd * v->dw * v->dh * sizeof(float));
		for (int j = 0; j < v->dh; j++)
		for (int i = 0; i < v->dw; i++)
		{
			// basically, redo the function "pan_repaint"
			static void (*win_to_img)(double p[2],struct pan_state*,double,double);
			win_to_img = window_to_image_apm;
			if (e->image_space)    win_to_img = window_to_image_apm;
			else if(e->force_exact) win_to_img = window_to_image_ex;
			int o = obtain_octave(e);
			int io = e->interpolation_order;
			//if (!e->image_space) continue;
			double p[2];
			win_to_img(p, e, i, j);
			if (p[0]<0 || p[1]<0 || p[0]>=v->w || p[1]>=v->h)
				continue;
			float pc = (p[0] + v->rgbiox) / 4;
			float qc = (p[1] + v->rgbioy) / 4;
			float g = i+0.1*j, c[pd];
			if(io==2)tiffo_getpixel_float_bilinear(&g,v->tg,0,*p,p[1]);
			else if(io==3)tiffo_getpixel_float_bicubic(&g,v->tg,0,*p,p[1]);
			else tiffo_getpixel_float_1     (&g,v->tg,0,p[0],p[1]);
			if (!v->gray_only) {
				if(io==2)tiffo_getpixel_float_bilinear(c,v->tc,0,pc,qc);
				else if(io==3)tiffo_getpixel_float_bicubic(c,v->tc,0,pc,qc);
				else tiffo_getpixel_float_1(c,v->tc,0,pc,qc);
			}
			buf_pan[i+j*v->dw] = g;
			for (int l = 0; l < pd; l++)
				buf_msi[(i+j*v->dw)*pd+l] = c[l];
		}
		iio_write_image_float_vec(fname_pan, buf_pan, v->dw, v->dh,  1);
		if (!v->gray_only)
			iio_write_image_float_vec(fname_msi, buf_msi, v->dw, v->dh, pd);
		fprintf(stderr, "dumped raw pan to file \"%s\"\n", fname_pan);
		if (!v->gray_only)
			fprintf(stderr, "dumped raw msi (pd=%d) to file \"%s\"\n",pd,fname_msi);
	}
	raw_collection_counter += 1;
}



// CALLBACK: pan_button_handler {{{1
static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	if (b == FTR_BUTTON_UP && m & FTR_MASK_CONTROL) {
		action_cycle_view(f, +1, x, y); return; }
	if (b == FTR_BUTTON_DOWN && m & FTR_MASK_CONTROL) {
		action_cycle_view(f, -1, x, y); return; }
	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
}


// CALLBACK: pan_motion_handler {{{1

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int unused_b, int m, int x, int y)
{
	(void)unused_b;
	struct pan_state *e = f->userdata;
	static double ox = 0, oy = 0;
	if (m & FTR_BUTTON_LEFT) action_offset_viewport(f, x - ox, y - oy);
	ox = x;
	oy = y;

	if (e->show_vertdir)
	{
		e->vdx = x;
		e->vdy = y;
		f->changed = 1;
	}
}

// CALLBACK: pan_key_handler {{{1
void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (k == '+' || k == '=') action_increase_zoom(f, f->w/2, f->h/2);
	if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);

	if (k == 'j') action_offset_rgbi(f, 0, -1);
	if (k == 'k') action_offset_rgbi(f, 0, 1);
	if (k == 'h') action_offset_rgbi(f, 1, 0);
	if (k == 'l') action_offset_rgbi(f, -1, 0);

	if (k == 'a') action_multiply_contrast(f, 1.3);
	if (k == 's') action_multiply_contrast(f, 1/1.3);
	if (k == 'u' && m & FTR_MASK_SHIFT){action_offset_base_h(f, +1);return;}
	if (k == 'd' && m & FTR_MASK_SHIFT){action_offset_base_h(f, -1);return;}
	if (k == 'u') {action_offset_base_h(f, +10);return;}
	if (k == 'd') {action_offset_base_h(f, -10);return;}

	if (k == 'e') action_toggle_exact_rpc(f);
	if (k == 'z') action_toggle_srtm4(f);
	if (k == 'v') action_toggle_srtm4_base(f);
	if (k == 'i') action_toggle_image_space(f, f->w/2, f->h/2);
	if (k == 'p') action_toggle_vertdir(f, x, y);
	if (k == '6') action_shift_so(f, -1);
	if (k == '7') action_shift_so(f, +1);
	if (k == 't') action_dump_view(f);
	if (k == '5') action_dump_collection(f);
	if (k == '4') action_dump_raw_fancy(f);
	if (k == '8') action_dump_raw_collection_fancy(f);
	if (k == ';') action_dump_raw(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// image flip operations
	if (k == ' ') action_cycle_view(f, 1, x, y);
	if (k == '\b') action_cycle_view(f, -1, x, y);
	if (k == '\'') action_toggle_diff_mode(f);
	if (k == 'y') action_toggle_log_scale(f);
	if (k == 'r') action_cycle_rotation_status(f);
	if (k == 'c') action_cycle_contrast(f, +1);
	if (k == 'C') action_cycle_contrast(f, -1);

	// arrows move the viewport
	if (k > 1000) {
		int d[2] = {0, 0};
		int inc = -BAD_MIN(f->w,f->h)/20;
		if (m & FTR_MASK_SHIFT  ) inc /= 10;
		if (m & FTR_MASK_CONTROL) inc *= 6;
		switch (k) {
		case FTR_KEY_LEFT : d[0] -= inc; break;
		case FTR_KEY_RIGHT: d[0] += inc; break;
		case FTR_KEY_UP   : d[1] -= inc; break;
		case FTR_KEY_DOWN : d[1] += inc; break;
		}
		if (k == FTR_KEY_PAGE_UP)   d[1] = +f->h/3;
		if (k == FTR_KEY_PAGE_DOWN) d[1] = -f->h/3;
		if (d[0] || d[1])
			action_offset_viewport(f, d[0], d[1]);
	}

	if (m & FTR_MASK_SHIFT) return;
	if (k=='1' ||k=='2' ||k=='3')
		action_select_interpolator(f, k-'0');
}

// CALLBACK: pan_resize {{{1
static void pan_resize(struct FTR *f, int k, int m, int x, int y)
{
	request_repaints(f);
}

// non-interactive control of the model {{{1
static int pan_non_interactive(struct pan_state *e, char *command_string)
{
	fprintf(stderr, "non-interactive rpcflip \"%s\"\n", command_string);

	// fill-in default options
	int x = e->view->w/2;
	int y = e->view->h/2;
	int w = 512;
	int h = 512;
	int base_h = 250;
	int interpord = 0;
	int qauto = 0;
	float zoom = 1;
	char outdir[FILENAME_MAX] = "/tmp/";
	char outnam[FILENAME_MAX] = "rpcflipo";

	// replace default options with the specified ones, if any
	char *delim = " \n\t,;", *tok = strtok(command_string, delim);
	do {
		fprintf(stderr, "\tmot = \"%s\"\n", tok);
		// change the fields that are recognized
		if (*tok == 'w') w = atoi(tok+1);
		if (*tok == 'h') h = atoi(tok+1);
		if (*tok == 'x') x = atoi(tok+1);
		if (*tok == 'y') y = atoi(tok+1);
		if (*tok == 'b') base_h = atoi(tok+1);
		if (*tok == 'i') interpord = atoi(tok+1);
		if (*tok == 'z') zoom = atof(tok+1);
		if (*tok == 'd') strncpy(outdir, tok+1, FILENAME_MAX);
		if (*tok == 'o') strncpy(outnam, tok+1, FILENAME_MAX);
		if (*tok == 'c') qauto = atoi(tok+1);
		if (*tok == 'p') {
			struct pan_view *v = e->view;
			if (!v->preview)v->preview = load_nice_preview(v->pfg,
							v->pfc, &v->pw, &v->ph);
			x = (atoi(tok+1) * e->view->w) / e->view->pw;
		}
		if (*tok == 'q') {
			struct pan_view *v = e->view;
			if (!v->preview) v->preview = load_nice_preview(v->pfg,
							v->pfc, &v->pw, &v->ph);
			y = (atoi(tok+1) * e->view->h) / e->view->ph;
		}
	} while ( (tok =strtok(NULL, delim)) );

	// dump output as requested
	struct FTR f[1]; // fake FTR, because the actions (wrongly) require it
	f->w = w; f->h = h; f->userdata = e;
	e->image_space = 1; // (this will be difficult to change)
	e->image_rotation_status = 0;
	e->interpolation_order = interpord;
	e->offset_x = x - e->image_space*w/(2*zoom);
	e->offset_y = y - e->image_space*h/(2*zoom);
	e->zoom_factor = zoom;
	e->base_h = base_h;
	e->qauto = qauto;
//#pragma omp parallel for // DOES NOT WORK (...leaving the domain of the sun)
	for (int i = 0; i < e->nviews; i++)
	{
		pan_repaint(e, w, h);
		struct pan_view *v = obtain_view(e);
		char buf[FILENAME_MAX];
		snprintf(buf, FILENAME_MAX, "%s/%s_%d.png", outdir, outnam, i);
		iio_write_image_uint8_vec(buf, v->display, w, h, 3);
		action_select_view(f, i+1, w/2, h/2); // not parallelizable
	}

	return 0;
}

// main {{{1
#include "pickopt.c"
int main_rpcflip_pan(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// process input arguments
	char *command_string = pick_option(&c, &v, "c", "");
	if (c < 5 || (c - 1) % 5 != 0) {
		fprintf(stderr, "usage:\n\t"
				"%s [P.TIF P.JPG P.RPC MS.TIF MS.JPG]+\n", *v);
		//                0  1     2     3     4      5
		return c;
	}
	int n = BAD_MIN((c - 1)/5,MAX_VIEWS);
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

	// if non-interactive, run the requested command string and exit
	if (*command_string)
		return pan_non_interactive(e, command_string);

	// open window
	//struct FTR f = ftr_new_window(320, 320);
	struct FTR f = ftr_new_window(800, 600);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , pan_key_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main_rpcflip_pank(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// process input arguments
	char *command_string = pick_option(&c, &v, "c", "");
	if (c < 4 || (c - 1) % 4 != 0) {
		fprintf(stderr, "usage:\n\t"
				"%s [P.TIF P.RPC MS.TIF MS.RPC]+\n", *v);
		//                0  1     2     3      4
		return c;
	}
	int n = BAD_MIN((c - 1)/4,MAX_VIEWS);
	fprintf(stderr, "we have %d views\n", n);
	char *filename_gtif[n];
	char *filename_grpc[n];
	char *filename_ctif[n];
	char *filename_crpc[n];
	for (int i = 0; i < n; i++)
	{
		filename_gtif[i] = v[4*i+1];
		filename_grpc[i] = v[4*i+2];
		filename_ctif[i] = v[4*i+3];
		filename_crpc[i] = v[4*i+4];
		fprintf(stderr, "view %d\n", i);
		fprintf(stderr, "\tgtif %s\n", filename_gtif[i]);
		fprintf(stderr, "\tgrpc %s\n", filename_grpc[i]);
		fprintf(stderr, "\tctif %s\n", filename_ctif[i]);
		fprintf(stderr, "\tcrpc %s\n", filename_crpc[i]);
	}

	// start state
	struct pan_state e[1];
	init_state_no_preview(e,
			filename_gtif,
			filename_grpc,
			filename_ctif,
			filename_crpc,
			n);

	// if non-interactive, run the requested command string and exit
	if (*command_string)
		return pan_non_interactive(e, command_string);

	// open window
	//struct FTR f = ftr_new_window(320, 320);
	struct FTR f = ftr_new_window(800, 600);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , pan_key_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main_rpcflip_gray(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// process input arguments
	char *command_string = pick_option(&c, &v, "c", "");
	if (c < 2 || (c - 1) % 2 != 0) {
		fprintf(stderr, "usage:\n\t"
				"%s [P.TIF P.RPC]+\n", *v);
		//                0  1     2
		return c;
	}
	int n = BAD_MIN((c - 1)/2,MAX_VIEWS);
	fprintf(stderr, "we have %d views\n", n);
	char *filename_gtif[n];
	char *filename_grpc[n];
	for (int i = 0; i < n; i++)
	{
		filename_gtif[i] = v[2*i+1];
		filename_grpc[i] = v[2*i+2];
		fprintf(stderr, "view %d\n", i);
		fprintf(stderr, "\tgtif %s\n", filename_gtif[i]);
		fprintf(stderr, "\tgrpc %s\n", filename_grpc[i]);
	}

	// start state
	struct pan_state e[1];
	init_state(e, filename_gtif, NULL, filename_grpc, NULL, NULL, n);

	// if non-interactive, run the requested command string and exit
	if (*command_string)
		return pan_non_interactive(e, command_string);

	// open window
	//struct FTR f = ftr_new_window(320, 320);
	struct FTR f = ftr_new_window(800, 600);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , pan_key_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main_rpcflip(int c, char *v[])
{
	if      (c > 1 && 0 == strcmp(v[1], "-g"))
		return main_rpcflip_gray(c - 1, v + 1);
	else if (c > 1 && 0 == strcmp(v[1], "-k"))
		return main_rpcflip_pank(c - 1, v + 1);
	else
		return main_rpcflip_pan(c, v);
}

int main(int c, char *v[])
{
	return main_rpcflip(c, v);
}

// vim:set foldmethod=marker:
