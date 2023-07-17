// CPU : the definitive interface for image processing practitioners

// cc -O2 cpu.c fancy_image.o iio.o ftr.o -o cpu -lX11 -lfftw3f -ltiff -ljpeg -lpng -lz -lm
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define CPU_SIGUSR
#ifdef CPU_SIGUSR
#include <signal.h>
#endif//CPU_SIGUSR

#include "fancy_image.h"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.h"

#define OMIT_MAIN_FONTU
#include "fontu.c" // todo: cherry-pick the required fontu functions
#include "fonts/xfonts_all.c"

#define OMIT_MAIN_SHADOWCAST
#include "shadowcast.c"

#define OMIT_MAIN_SARSIM
#include "sarsim.c"

#include "xmalloc.c"

#define WHEEL_FACTOR 2
#define MAX_PYRAMID_LEVELS 30

#define MAX_IMAGES 1000

// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data (for the current image)
	int w, h;
	struct fancy_image *i;

	// 2. view port parameters
	double zoom_factor, offset_x, offset_y;
	double a, b, bbb[3];

	// 3. ancillary data, viewing configuration
	float chan_mu[3];
	float chan_size[3];
	float p;
	bool auto_qauto;

	// 4. roi
	int roi; // 0=nothing 1=dftwindow 2=rawfourier 3=ppsmooth 4=circppsm
	int roi_x, roi_y, roi_w;

	// 5. user interface
	struct bitmap_font font[5];
	int hud;

	// 6. topographic mode
	int topographic_mode;
	// 0=no, 1=botw, 2=shadow, 3=linear, 4=lambert, 5=specular
	float topographic_sun[3];
	float topographic_scale;
	float topographic_P;
	float topographic_spread;

	// 7. actual image data for the whole series
	int i_idx, i_num;
	struct fancy_image *i_tab[MAX_IMAGES];
	char *i_name[MAX_IMAGES];
};

// change of coordinates: from window "int" pixels to image "double" point
static void window_to_image(double p[2], struct pan_state *e, int i, int j)
{
	p[0] = e->offset_x + i / e->zoom_factor;
	p[1] = e->offset_y + j / e->zoom_factor;
}

// change of coordinates: from image "double" point to window "int" pixel
static void image_to_window(int i[2], struct pan_state *e, double x, double y)
{
	i[0] = floor((x - e->offset_x) * e->zoom_factor );
	i[1] = floor((y - e->offset_y) * e->zoom_factor );
}

// TODO: refactor this function with "pixel" and "colormap3"
// add various shader options
static void get_rgb_from_vec(float *rgb, struct pan_state *e, float *vec)
{
	rgb[0] = rgb[1] = rgb[2] = vec[0];
	if (e->i->pd > 1)
		rgb[1] = rgb[2] = vec[1];
	if (e->i->pd > 2)
		rgb[2] = vec[2];
	if (e->i->pd == 4)
		rgb[1] = rgb[1] + 0.1 * vec[3];
}

//static float fancy_interpolate(struct fancy_image *f, int oct,
//		float p, float q, int l
//{
//}

// like n%p, but works for all numbers
static int gmod(int x, int m)
{
	int r = x % m;
	return r < 0 ? r + m : r;
}

// evaluate the value a position (p,q) in image coordinates
static void pixel(float *out, struct pan_state *e, double p, double q)
{
	//if(p<0||q<0){out[0]=out[1]=out[2]=170;return;}// TODO: kill this
	//if(p>=e->i->w||q>=e->i->h){out[0]=out[1]=out[2]=85;return;}
	if (p < 0 || q < 0 || p >= e->i->w || q >= e->i->h) {
		int ip = p+256;
		int iq = q+256;
		int pip = gmod(ip/256, 2);
		int piq = gmod(iq/256, 2);
		int val = gmod(pip+piq,2);
		out[0] = out[1] = out[2] = 127+val*64;
		return;
	}

	int oct = 0;
	if (e->zoom_factor < 0.9999)
	{
		int s = round(log2(1/e->zoom_factor));
		if (s < 0) s = 0;
		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
		int sfac = 1<<(s);
		oct = s;
		p /= sfac;
		q /= sfac;
	}
	for (int i = 0; i < e->i->pd; i++)
	{
		// TODO, interpolation
		out[i] = fancy_image_getsample_oct(e->i, oct, p, q, i);
	}
}

// evaluate the value a position (p,q) in image coordinates
static float pixel_height(struct pan_state *e, double p, double q)
{
	if (p < 0 || q < 0 || p >= e->i->w || q >= e->i->h)
		return NAN;

	int oct = 0;
	if (e->zoom_factor < 0.9999)
	{
		int s = round(log2(1/e->zoom_factor));
		if (s < 0) s = 0;
		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
		int sfac = 1<<(s);
		oct = s;
		p /= sfac;
		q /= sfac;
	}

	return fancy_image_getsample_oct(e->i, oct, p, q, 0);
}

static void pixel_rgbf(float out[3], struct pan_state *e, double p, double q)
{
	float v[e->i->pd];
	v[e->i->pd-1] = 0; // remove an idiotic gcc warning
	pixel(v, e, p, q);
	get_rgb_from_vec(out, e, v);
}

static void action_print_value_under_cursor(struct FTR *f, int x, int y)
{
	if (x<f->w && x>=0 && y<f->h && y>=0) {
		struct pan_state *e = f->userdata;
		double p[2];
		window_to_image(p, e, x, y);
		float c[3];
		pixel_rgbf(c, e, p[0], p[1]);
		//interpolate_at(c, e->frgb, e->w, e->h, p[0], p[1]);
		//printf("%g\t%g\t: %g\t%g\t%g\n", p[0], p[1], c[0], c[1], c[2]);
		printf("not implemented %g %g : %g %g %g\n",
				p[0], p[1], c[0], c[1], c[2]);
	}
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->zoom_factor = 1;
	e->offset_x = 0;
	e->offset_y = 0;
	e->a = 1;
	e->b = 0;
	e->bbb[0] = e->bbb[1] = e->bbb[2] = 0;

	e->p = INFINITY;
	e->auto_qauto = false;
	e->hud = 0;

	e->roi = 0;
	e->roi_x = f->w / 2;
	e->roi_y = f->h / 2;
	e->roi_w = 73; // must be odd

	e->topographic_mode = 0;
	e->topographic_sun[0] = 1/sqrt(3);
	e->topographic_sun[1] = 1/sqrt(3);
	e->topographic_sun[2] = 1/sqrt(3);
	e->topographic_scale = 1;
	e->topographic_P = 30;

	f->changed = 1;
}

static void action_contrast_change(struct FTR *f, float afac, float bshift)
{
	struct pan_state *e = f->userdata;

	e->a *= afac;
	e->b += bshift;

	f->changed = 1;
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

#include "smapa.h"
SMART_PARAMETER(SIGMAFACTOR_P1,2)
SMART_PARAMETER(SIGMAFACTOR_P2,3)

void compute_scalar_position_and_size(
		float *out_pos, float *out_siz,
		float *x, int n,
		float p)
{
	fprintf(stderr, "ska_pos_s    p=%g     len=%d    first=%g\n",
			p, n, x[0]);

	if (false) ;
	else if (isnan(p))
	{
		*out_pos = 127.5;
		*out_siz = 255;
	}
	else if (p == INFINITY)
	{
		float min = INFINITY;
		float max = -INFINITY;
		for (int i = 0; i < n; i++)
		{
			min = fmin(min, x[i]);
			max = fmax(max, x[i]);
		}
		*out_pos = (min + max) / 2;
		*out_siz = max - min;
	}
	else if (p == 2) // mean, average
	{
		long double mu = 0;
		long double sigma = 0;
		for (int i = 0; i < n; i++)
			mu += x[i];
		mu /= n;
		for (int i = 0; i < n; i++)
			sigma = hypot(sigma, x[i] - mu);
		sigma /= sqrt(n);
		*out_pos = mu;
		*out_siz = sigma * SIGMAFACTOR_P2();
	}
	else if (p == 1)
	{
		float *t = xmalloc(n * sizeof*t);
		for (int i = 0; i < n; i++)
			t[i] = x[i];
		qsort(t, n, sizeof*t, compare_floats);
		float med = t[n/2];
		free(t);
		long double aad = 0;
		for (int i = 0; i < n; i++)
			aad += fabs(x[i] - med);
		aad /= n;
		*out_pos = med;
		*out_siz = aad * SIGMAFACTOR_P1();
	}
	else fail("unrecognized p=%g\n", p);
	// todo: implement other measures

	fprintf(stderr, "\tpos size = %g %g\n", *out_pos, *out_siz);
}

void compute_split_position_and_size(
		float *m, // array of output mus
		float *s, // array of output sigmas
		float *x, // input array of n d-dimsensional vectors
		int n,    // number of input vectors
		int d,    // dimension of each vector
		float p   // p robustness parameter (p=1,2,inf)
		)
{
	// we do the coward thing, component by component
	float *t = xmalloc(n * sizeof *t);
	for (int l = 0; l < d; l++)
	{
		int N = 0;
		for (int i = 0; i < n; i++)
			if (isfinite(x[i*d+l]))
				t[N++] = x[i*d+l];
		compute_scalar_position_and_size(m + l, s + l, t, N, p);
	}
	free(t);
}

static void action_qauto(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	float m = INFINITY, M = -m;
	int pid = 3;
	m = 0;
	M = 255;
	//for (int i = 0; i < 3 * e->pyr_w[pid] * e->pyr_h[pid]; i++)
	//{
	//	float g = e->pyr_rgb[pid][i];
	//	m = fmin(m, g);
	//	M = fmax(M, g);
	//}

	e->a = 255 / ( M - m );
	e->b = 255 * m / ( m - M );
	e->bbb[0] = e->b;
	e->bbb[1] = e->b;
	e->bbb[2] = e->b;

	f->changed = 1;
}

static void action_qauto2(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	//static int qauto_p_idx = 0;
	//float p[4] = {NAN, INFINITY, 2, 1};
	//qauto_p_idx = (qauto_p_idx + 1) % 4;
	//fprintf(stderr, "qauto p = %g\n", p[qauto_p_idx]);

	// build an array of pixel values from the current screen
	// TODO: cache this shit, we are computing this twice!
	int pd = e->i->pd;
	float *t = xmalloc(f->w * f->h * pd * sizeof*t);
	int winsize = 0;
	for (int j = 0.1*f->h; j < 0.9*f->h; j++)
	for (int i = 0.1*f->w; i < 0.9*f->w; i++)
	{
		double p[2]; window_to_image(p, e, i, j);
		pixel(t + pd * winsize++, e, p[0], p[1]);
	}

	// extract mu/sigma for each channel
	float mu[pd], sigma[pd];
	compute_split_position_and_size(mu, sigma, t, winsize, pd, e->p);
	free(t);

	for (int i = 0; i < pd; i++)
		fprintf(stderr, "musigma[%d] = %g %g\n", i, mu[i], sigma[i]);


	// adapt mu/sigma to color
	float mu_rgb[3], sigma_rgb[3];
	get_rgb_from_vec(mu_rgb   , e, mu   );
	get_rgb_from_vec(sigma_rgb, e, sigma);

	// change contrast viewport accordingly
	e->a = 255 * 3 / (sigma_rgb[0] + sigma_rgb[1] + sigma_rgb[2]);
	e->bbb[0] = 127.5 - e->a * mu_rgb[0];
	e->bbb[1] = 127.5 - e->a * mu_rgb[1];
	e->bbb[2] = 127.5 - e->a * mu_rgb[2];

	f->changed = 1;
}

static void action_toggle_roi(struct FTR *f, int x, int y, int dir)
{
	struct pan_state *e = f->userdata;
	e->roi = (e->roi + (dir?-1:1)) % 3;
	fprintf(stderr, "ROI SWITCH(%d) = %d\n", dir, e->roi);
	e->roi_x = x - e->roi_w / 2;
	e->roi_y = y - e->roi_w / 2;
	f->changed = 1;
}

static void action_cycle_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->hud = !e->hud;
	f->changed = 1;
}

static void action_toggle_topography(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->topographic_mode = (1 + e->topographic_mode) % 7;
	f->changed = 1;
}

static void action_topography_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;
	e->topographic_scale *= factor;
	fprintf(stderr, "TOPOGRAPHIC SCALE = %g\n", e->topographic_scale);
	f->changed = 1;
}

static void action_topography_Pspan(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;
	e->topographic_P *= factor;
	fprintf(stderr, "TOPOGRAPHIC P = %g\n", e->topographic_P);
	f->changed = 1;
}


static void action_toggle_aqauto(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->auto_qauto = !e->auto_qauto;
	if (e->auto_qauto)
		action_qauto2(f);
}

static void action_toggle_p(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	if (false) ;
	else if (isnan(e->p)) e->p = 1;
	else if (e->p == 1) e->p = 2;
	else if (e->p == 2) e->p = INFINITY;
	else if (e->p == INFINITY) e->p = NAN;

	fprintf(stderr, "P = %g\n", e->p);
	action_qauto2(f);
	f->changed = 1;
}

static void action_roi_embiggen(struct FTR *f, float s)
{
	struct pan_state *e = f->userdata;
	e->roi_w *= s;
	fprintf(stderr, "ROI %d\n", e->roi_w);
	f->changed = 1;
}


static void action_move_roi(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;
	if (e->roi)
	{
		e->roi_x = x - e->roi_w / 2;
		e->roi_y = y - e->roi_w / 2;
		f->changed = 1;
	}
}

static void action_move_sun(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;
	float R = fmin(f->w, f->h)/2;
	float p = (x - f->w/2) / R;
	float q = (y - f->h/2) / R;
	float r = sqrt(1-p*p-q*q);
	if (!isfinite(r)) r = 0;
	float n = hypot(p, hypot(q, r));
	e->topographic_sun[0] = p/n;
	e->topographic_sun[1] = q/n;
	e->topographic_sun[2] = r/n;
	fprintf(stderr, "SUN = %g %g %g\n", p/n, q/n, r/n);
	f->changed = 1;
}

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	float c[3];
	pixel_rgbf(c, e, p[0], p[1]);
	float C = (c[0] + c[1] + c[2])/3;

	e->bbb[0] = 127.5 - e->a * c[0];
	e->bbb[1] = 127.5 - e->a * c[1];
	e->bbb[2] = 127.5 - e->a * c[2];

	e->b = 127.5 - e->a * C;

	f->changed = 1;
}

static void action_contrast_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->b)/ e->a;
	float ccc[3];
	for(int l=0;l<3;l++) ccc[l] = (127.5 - e->bbb[l]) / e->a;
	e->a *= factor;
	e->b = 127.5 - e->a * c;
	for(int l=0;l<3;l++) e->bbb[l] = 127.5 - e->a * ccc[l];

	f->changed = 1;
}

static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double c[2];
	window_to_image(c, e, x, y);

	e->zoom_factor *= F;
	e->offset_x = c[0] - x/e->zoom_factor;
	e->offset_y = c[1] - y/e->zoom_factor;
	fprintf(stderr, "\t zoom changed %g\n", e->zoom_factor);

	if (e->auto_qauto) action_qauto2(f);
	f->changed = 1;
}

static void action_reset_zoom_only(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	action_change_zoom_by_factor(f, x, y, 1/e->zoom_factor);
}

static void action_increase_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, WHEEL_FACTOR);
}

static void action_decrease_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, 1.0/WHEEL_FACTOR);
}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->zoom_factor;
	e->offset_y -= dy/e->zoom_factor;

	if (e->auto_qauto) action_qauto2(f);

	f->changed = 1;
}

static void action_reload_image(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	fancy_image_reload(e->i);
	f->changed = 1;
}

#ifdef CPU_SIGUSR
static struct FTR *global_f_for_sigusr = NULL;
void handle_signal(int n)
{
	void ftr_x11_force_redraw(struct FTR*);
	if (false) {
		;
	} else if (n == SIGUSR1) {
		fprintf(stderr, "signal %d received\n", n);
		action_reload_image(global_f_for_sigusr);
#if FTR_BACKEND == 'x'
		ftr_x11_force_redraw(global_f_for_sigusr);
#endif
	} else if (n == SIGUSR2) {
		fprintf(stderr, "signal %d received\n", n);
		ftr_notify_the_desire_to_stop_this_loop(global_f_for_sigusr, 1);
#if FTR_BACKEND == 'x'
		ftr_x11_force_redraw(global_f_for_sigusr);
#endif
	}
	else
		fprintf(stderr, "signal %d unrecognized\n", n);

}
#endif//CPU_SIGUSR


static void action_update_window_title(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	int n = 2*FILENAME_MAX;
	char t[n];
	if (e->i_num == 1)
		snprintf(t, n, "%s", e->i_name[0]);
	else
		snprintf(t, n, "%s [%d/%d]", e->i_name[e->i_idx],
				1 + e->i_idx, e->i_num);
	ftr_change_title(f, t);
}

static void action_flipn(struct FTR *f, int idx)
{
	struct pan_state *e = f->userdata;
	if (idx < 0 || idx >= e->i_num || idx == e->i_idx)
		fprintf(stderr, "warning: no image to flip\n");
	else {
		fprintf(stderr, "flip %d to %d\n", e->i_idx, idx);
		e->i = e->i_tab[idx];
		e->i_idx = idx;
		e->w = e->i->w;
		e->h = e->i->h;
		f->changed = 1;
	}
}

static void action_flip(struct FTR *f, int o)
{
	struct pan_state *e = f->userdata;
	int i = e->i_idx;
	i = i + o;
	if (i < 0) i = e->i_num - 1;
	if (i >= e->i_num) i = 0;
	action_flipn(f, i);
	action_update_window_title(f);
}

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_cpu_%d.png", c);
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(n, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}

static bool insideP(int w, int h, int i, int j)
{
	return i>=0 && j>=0 && i<w && j<h;
}

#define OMIT_BLUR_MAIN
#include "blur.c"

#define OMIT_PPSMOOTH_MAIN
#include "ppsmooth.c"

static void ppsmooth_vec(float *y, float *x, int w, int h, int pd)
{
	float *x_split = xmalloc(w*h*pd*sizeof*x);
	float *y_split = xmalloc(w*h*pd*sizeof*x);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		x_split[l*w*h+(j*w+i)] = x[pd*(j*w+i)+l];

	ppsmooth_split(y_split, x_split, w, h, pd);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		y[pd*(j*w+i)+l] = y_split[l*w*h+(j*w+i)];

	free(x_split);
	free(y_split);
}

// note: assumes 3-dimensional pixels
static void transform_roi_buffers_old(float *y, float *x, int n)
{
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
		float *xx = x + 3*(j*n + i);
		float *yy = y + 3*(j*n + i);
		float g = hypot(xx[0], hypot(xx[1],xx[2]))/sqrt(3);
		yy[0] = sqrt(g/255)*255;
		yy[1] = g/2;
		yy[2] = g;
	}
	if (false) return;
	float x_p[3*n*n];
	float x_s[3*n*n];
	ppsmooth_vec(x_p, x, n, n, 3);
	for (int i = 0; i < 3*n*n; i++)
		x_s[i] = x[i] - x_p[i];
	float param[1] = {4};
	blur_2d(y, x_p, n, n, 3, "laplace", param, 1);
	for (int i = 0; i < 3*n*n; i++)
		y[i] += x_s[i];
}

static void circular_roi(float *y, float *x, int n)
{
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
	{
		if (hypot(i-0.5*n+0.5,j-0.5*n+0.5) < 0.5*n-4)
			y[(j*n+i)*3+l] = NAN;
		else
			y[(j*n+i)*3+l] = x[(j*n+i)*3+l];
	}
	simplest_inpainting_vec(y, n, n, 3);
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	for (int l = 0; l < 3; l++)
		if (hypot(i-0.5*n+0.5,j-0.5*n+0.5) < 0.5*n-4)
			y[(j*n+i)*3+l] = x[(j*n+i)*3+l] - y[(j*n+i)*3+l];// + 127;
}

static void fourier_roi(float *y, float *x, int n)
{
	float *c = xmalloc(n*n*sizeof*c);
	float *ys = xmalloc(n*n*sizeof*c);
	fftwf_complex *fc = fftwf_xmalloc(n*n*sizeof*fc);
	for (int l = 0; l < 3; l++)
	{
		for (int i = 0; i < n*n; i++)
			c[i] = x[3*i+l];
		fft_2dfloat(fc, c, n, n);
		for (int i = 0; i < n*n; i++)
			//ys[i] = cabs(fc[i]);
			ys[i] = 255*(log(cabs(fc[i])/255)+0.5)/5;
		//float fac = n*13;
		float fac = n*3;
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			int ii = (i + n/2) % n;
			int jj = (j + n/2) % n;
			float norm = hypot(i-n/2-1, j-n/2-1) / fac;
			y[3*(j*n+i)+l] = ys[jj*n+ii] * 1;//norm;
		}
	}
	free(ys);
	free(c);
	fftwf_free(fc);
	//float param[1] = {0.5};
	//blur_2d(y, y, n, n, 3, "cauchy", param, 1);
}

static void autocorrelation_roi(float *y, float *x, int n)
{
	float *c = xmalloc(n*n*sizeof*c);
	float *ys = xmalloc(n*n*sizeof*c);
	fftwf_complex *fc = fftwf_xmalloc(n*n*sizeof*fc);
	for (int l = 0; l < 3; l++)
	{
		for (int i = 0; i < n*n; i++)
			c[i] = x[3*i+l];
		fft_2dfloat(fc, c, n, n);
		for (int i = 0; i < n*n; i++)
			fc[i] = cabs(fc[i]);
		ifft_2dfloat(c, fc, n, n);
		for (int i = 0; i < n*n; i++)
			ys[i] = c[i];
		//float fac = n*13;
		float fac = n*3;
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			int ii = (i + n/2) % n;
			int jj = (j + n/2) % n;
			//float norm = hypot(i-n/2-1, j-n/2-1) / fac;
			y[3*(j*n+i)+l] = ys[jj*n+ii] * 1;//norm;
		}
	}
	free(ys);
	free(c);
	fftwf_free(fc);
}

static void transform_roi_buffers(float *y, float *x, int n, int roi)
{
	//if (roi == 4) { circular_roi(y, x, n); return; }
	float x_p0[3*n*n], *x_p = x_p0;
	ppsmooth_vec(x_p0, x, n, n, 3);
	if (roi == 1) { fourier_roi(y, x_p, n); return; }
	if (roi == 2) { autocorrelation_roi(y, x_p, n); return; }
	//if (roi == 2) x_p = x;
	//if (roi == 3) { for (int i = 0; i < 3*n*n; i++) y[i]=x_p0[i]; return; }
}

// TODO: combine "colormap3" and "pixel"
static void colormap3(unsigned char *rgb, struct pan_state *e, float *frgb)
{
	for (int l = 0; l < 3; l++)
	{
		if (!isfinite(frgb[l]))
			//rgb[l] = l==1?0:255; // show NAN as magenta
			rgb[l] = l==2?127:0; // show NAN as dark blue
		else {
			//float g = e->a * c[l] + e->b;
			float g = e->a * frgb[l] + e->bbb[l];
			if      (g < 0)   rgb[l] = 0  ;
			else if (g > 255) rgb[l] = 255;
			else              rgb[l] = g  ;
		}
	}
}

static struct bitmap_font *get_font_for_zoom(struct pan_state *e, double z)
{
	return e->font + 1 * (z > 50) + 3 * (z > 100);
}

static void expose_pixel_values(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	double zf = e->zoom_factor;
	if (zf > 30) // "zoom_factor equals the pixel size"
	{
		// find first inner pixel in the image domain
		double p[2];
		window_to_image(p, e, 0, 0);
		double ip[2] = {ceil(p[0])+0.5, ceil(p[1])+0.5}; // port coords
		int ij[2]; // window coords
		image_to_window(ij, e, ip[0], ip[1]);
		for (int jj = ij[1] - zf; jj < f->h + zf; jj += zf)
		for (int ii = ij[0] - zf; ii < f->w + zf; ii += zf)
		{
			window_to_image(p, e, ii, jj);
			float c[3]; pixel_rgbf(c, e, p[0], p[1]);
			uint8_t rgb[3]; colormap3(rgb, e, c);
			uint8_t color_white[3] = {255, 255, 255};
			uint8_t color_green[3] = {0, 255, 0};
			uint8_t color_red[3] = {255, 0, 0};
			uint8_t bg[3] = {0, 0, 0};
			uint8_t *fg = rgb[1] < rgb[0] ? color_green : color_red;
			if (hypot(rgb[0], rgb[1]) > 250 && rgb[0] < 2*rgb[1])
				fg = color_red;
			char buf[300];
			int l = 0;
			if (e->i->pd > 2 && zf > 50)
			//if (e->zoom_factor>100 && (c[0]!=c[1] || c[0]!=c[2]))
				l=snprintf(buf,300,"%g\n%g\n%g",c[0],c[1],c[2]);
			else if (e->i->pd == 2 && zf > 50)
				l=snprintf(buf,300,"%g\n%g",c[0],c[1]);
			else if (e->i->pd == 1)
			//if (c[0] == c[1] && c[0] == c[2])
				l=snprintf(buf, 300, "%g", c[0]);
			struct bitmap_font *font = get_font_for_zoom(e, zf);
			uint8_t *rbg = zf>100?bg:NULL;
			if (l) put_string_in_rgb_image(f->rgb, f->w, f->h,
					ii-2.5*font->width, jj-1.5*font->height,
					fg, rbg, 0, font, buf);
			if (zf>129)
			{
				snprintf(buf,300,"%g %g\n",
						floor(p[0]), floor(p[1]));
				put_string_in_rgb_image(f->rgb, f->w, f->h,
					ii - zf/2,
					jj - zf/2,
					color_white,rbg,0,font,buf);
			}
		}

	}
}

static void expose_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	uint8_t fg[3] = {0, 0, 0};
	uint8_t bg[3] = {255, 255, 255};
	put_string_in_rgb_image(f->rgb, f->w, f->h, 10, 10, fg, bg, 0,
			e->font+4, e->i_name[e->i_idx]);
}

static void expose_roi(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	float buf_in [3 * e->roi_w * e->roi_w];
	float buf_out[3 * e->roi_w * e->roi_w];
	for (int j = 0; j < e->roi_w; j++)
	for (int i = 0; i < e->roi_w; i++)
	{
		// i,h       = position inside ROI
		// ii,jj     = position inside viewport
		// p[0],p[1] = position inside image domain
		int ii = i + e->roi_x - e->roi_w/2;
		int jj = j + e->roi_y - e->roi_w/2;
		double p[2];
		window_to_image(p, e, ii, jj);
		float *c = buf_in + 3 * (j * e->roi_w + i);
		pixel_rgbf(c, e, p[0], p[1]);
	}
	transform_roi_buffers(buf_out, buf_in, e->roi_w, e->roi);
	float A = e->a;
	float BBB[3] = {e->bbb[0], e->bbb[1], e->bbb[2]};
	if (4 == e->roi) {
		BBB[0] = BBB[1] = BBB[2] = 127;
	}
	for (int j = 0; j < e->roi_w; j++)
	for (int i = 0; i < e->roi_w; i++)
	{
		int ii = i + e->roi_x - e->roi_w/2;
		int jj = j + e->roi_y - e->roi_w/2;
		float *c = buf_out + 3 * (j * e->roi_w + i);
		if (insideP(f->w, f->h, ii, jj))
		{
			if (4 == e->roi &&
				hypot(i-0.5*e->roi_w+0.5,
					j-0.5*e->roi_w+0.5) >= 0.5*e->roi_w-4)
				continue;
			uint8_t *cc = f->rgb + 3 * (jj * f->w + ii);
			for (int l = 0; l < 3; l++)
			{
				float g = A * c[l] + BBB[l];
				if      (g < 0)   cc[l] = 0  ;
				else if (g > 255) cc[l] = 255;
				else              cc[l] = g  ;
				//cc[l] = c[l];
			}
		}
	}
}

static unsigned char bclamp(float x)
{
	if (x < 0) return 0;
	if (x > 255) return 255;
	return x;
}

static void getpercentiles(float *m, float *M, float *x, int n, float q)
{
	float *t = xmalloc(n * sizeof*t);
	memcpy(t, x, n*sizeof*t);
	int N = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(x[i]))
			t[N++] = x[i];
	qsort(t, N, sizeof*t, compare_floats);
	int a = q/100*N;
	int b = (1-q/100)*N;
	if (a < 0) a = 0;
	if (b < 0) b = 0;
	if (a >= N) a = N-1;
	if (b >= N) b = N-1;
	*m = t[a];
	*M = t[b];
	free(t);
}

// rewrites x a bit (e.g. to fill-in nans)
static void colorize_botw(uint8_t *y, float *x, int w, int h)
{
	// apply botw palette
	// TODO: add some 2% saturation here
	float m, M;
	getpercentiles(&m, &M, x, w*h, 1.0);
	for (int i = 0; i < w*h; i++)
	if (isfinite(x[i]) && x[i] != 32768) {
		float t = (x[i] - m) / (M - m);
		y[3*i+0] = bclamp( (1 - t)*65 + t*200 );
		y[3*i+1] = bclamp( (1 - t)*49 + t*200 );
		y[3*i+2] = bclamp( (1 - t)*5  + t*180 );
	} else {
		y[3*i+0] = 20;
		y[3*i+1] = 100;
		y[3*i+2] = 255;
		x[i] = 0;
	}

	// ppsmooth
	float *s = xmalloc(w*h*sizeof*s);
	ppsmooth(s, x, w, h);

	// compute lssao shading
	float *z = xmalloc(w * h * sizeof*z);
	fftwf_complex *S = fftwf_xmalloc(w * h * sizeof*S);
	fft_2dfloat(S, s, w, h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = i < w/2 ? i : i - w;
		int jj = j < h/2 ? j : j - h;
		S[j*w+i] *= hypot(ii*h, jj*w);
	}
	ifft_2dfloat(z, S, w, h);
	fftwf_free(S);
	for (int i = 0; i < w*h; i++)
		if (z[i] > 0)
			z[i] /= 4;
	// combine shading and palette
	getpercentiles(&m, &M, z, w*h, 5.0);
	for (int i = 0; i < w*h; i++)
		z[i] =  (z[i] - m) / (M - m);
	for (int i = 0; i < w*h; i++)
	for (int k = 0; k < 3; k++)
		y[3*i+k] = bclamp(y[3*i+k] * (0.5+z[i]/2));
		//y[3*i+k] = bclamp(255*z[i]);
	free(z);
}

static void expose_topography(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	if (e->topographic_mode == 1) // botw
	{
		float   *x = xmalloc(1 * f->w * f->h * sizeof*x);
		uint8_t *y = xmalloc(3 * f->w * f->h * sizeof*y);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		{
			double p[2];
			window_to_image(p, e, i, j);
			x[j*f->w+i] = pixel_height(e, p[0], p[1]
					) / e->topographic_scale;
		}
		colorize_botw(y, x, f->w, f->h);
		for (int i = 0; i < 3 * f->w * f->h; i++)
			f->rgb[i] = y[i];
		free(x);
		free(y);
		return;
	}

	if (e->zoom_factor != 1 || e->i->pd != 1) return;

	if (e->topographic_mode == 2) // shadows
	{
		float *x = xmalloc(f->w * f->h * sizeof*x);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
			x[j*f->w+i] = fancy_image_getsample(
					e->i,
					i + e->offset_x,
					j + e->offset_y,
					0) / e->topographic_scale;
		cast_shadows(x, f->w, f->h,
				-e->topographic_sun[0],
				-e->topographic_sun[1],
				tan(asin(e->topographic_sun[2])));
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		for (int k = 0; k < 3; k++)
			f->rgb[(j*f->w+i)*3+k] = 255*isfinite(x[f->w*j+i]);
		free(x);
		return;
	}


	if (e->topographic_mode == 6) // radar
	{
		float *x = xmalloc(f->w * f->h * sizeof*x);
		float *y = xmalloc(f->w * f->h * sizeof*x);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
			x[j*f->w+i] = fancy_image_getsample(
					e->i,
					i + e->offset_x,
					j + e->offset_y,
					0) / e->topographic_scale;
		float a = -acos(e->topographic_sun[0]);
		if (e->topographic_sun[0] < 0)
			a += M_PI;
		radar_sim_horizontal(y, x, f->w, f->h, a);
		//fprintf(stderr, "y[0,200] = %g\n", y[(f->w*200)+0]);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		for (int k = 0; k < 3; k++)
			//f->rgb[(j*f->w+i)*3+k] = e->topographic_P*y[f->w*j+i];
			f->rgb[(j*f->w+i)*3+k] = bclamp(e->topographic_P*y[f->w*j+i]);
		//{
		//	f->rgb[(j*f->w+i)*3+0] = e->topographic_P*y[f->w*j+i];
		//	f->rgb[(j*f->w+i)*3+1] = e->topographic_P*y[f->w*j+i];
		//	f->rgb[(j*f->w+i)*3+2] = e->topographic_P*y[f->w*j+i];
		//}
		//fprintf(stderr, "frgb[0,200] = %d %d %d\n",
		//		f->rgb[(f->w*200)*3+0],
		//		f->rgb[(f->w*200)*3+1],
		//		f->rgb[(f->w*200)*3+2]);
		free(x);
		free(y);
		return;
	}

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		int p = e->offset_x + i;
		int q = e->offset_y + j;
		float s = e->topographic_scale;
		float h = fancy_image_getsample(e->i, p, q, 0) / s;
		float h10 = fancy_image_getsample(e->i, p+1, q, 0) / s;
		float h01 = fancy_image_getsample(e->i, p, q+1, 0) / s;
		float hx = (h10 - h);
		float hy = (h01 - h);
		float *S = e->topographic_sun;

		float c = 0;
		unsigned char *rgb = f->rgb + 3 * (j * f->w + i);
		switch(e->topographic_mode) {
		//case 1: // shadows
			//break;
		case 3: // linearized lambertian
			c = - hx * S[0] - hy * S[1];
			rgb[0] = rgb[1] = rgb[2] = bclamp(127 + 40 * c);
			break;
		case 4: // lambertian (Gouraud)
			c = S[2] - hx * S[0] - hy * S[1];
			c /= sqrt(1 + hx*hx + hy*hy);
			rgb[0] = rgb[1] = rgb[2] = bclamp(e->topographic_P
					* fmax(0, c));
			break;
		case 5: { // specular (Blinn-Phong)
			float N[3] = {-hx, -hy, 1};
			float n = hypot(N[2], hypot(N[1], N[0]));
			N[0]/=n; N[1]/=n; N[2]/=n; // N = unit normal
			float H[3] = {S[0], S[1], S[2]+1};
			n = hypot(H[2], hypot(H[1], H[0]));
			H[0]/=n; H[1]/=n; H[2]/=n; // H = half-angle direction
			float k = pow(
					fmax(0, H[0]*N[0]+H[1]*N[1]+H[2]*N[2]),
					e->topographic_P);
			rgb[0] = rgb[1] = rgb[2] = bclamp(255 * k);
			break; }
		//case 6: // radar
			//break;
		default: fail("impossible topographic condition %d",
					 e->topographic_mode);
		}
	}
}

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->topographic_mode)
	{
		expose_topography(f);
		f->changed = 1;
		goto cont;
	}

	// expose the whole image
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		float c[3];
		pixel_rgbf(c, e, p[0], p[1]);
		unsigned char *cc = f->rgb + 3 * (j * f->w + i);
		colormap3(cc, e, c);
	}

cont:
	// if pixels are "huge", show their values
	if (e->zoom_factor > 30)
		expose_pixel_values(f);

	// if HUD, expose the hud
	if (e->hud)
		expose_hud(f);

	// if ROI, expose the roi
	if (e->roi && !e->topographic_mode)
		expose_roi(f);

	// mark shit as changed
	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "motion b=%d m=%d (%d %d)\n", b, m, x, y);
	struct pan_state *e = f->userdata;

	static double ox = 0, oy = 0;

	if (m & FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	if (m & FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (m & FTR_MASK_SHIFT)    action_center_contrast_at_point(f, x, y);
	if (m & FTR_MASK_SHIFT && e->topographic_mode) action_move_sun(f, x, y);

	action_move_roi(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "button b=%d m=%d\n", b, m);
	struct pan_state *e = f->userdata;

	if (e->roi && b == FTR_BUTTON_UP) { action_roi_embiggen(f,1.1); return; }
	if (e->roi && b == FTR_BUTTON_DOWN){action_roi_embiggen(f,1/1.1); return; }
	if (b == FTR_BUTTON_UP && m & FTR_MASK_SHIFT) {
		action_contrast_span(f, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && m & FTR_MASK_SHIFT) {
		action_contrast_span(f, 1.3); return; }
	if (b == FTR_BUTTON_RIGHT && m & FTR_MASK_CONTROL) {
		action_reset_zoom_only(f, x, y); return; }
	if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
}

static void key_handler_print(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "key pressed %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
}

static void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);
	//fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
	//		k, isprint(k)?k:' ', m, x, y);

	//if (k == '+') action_increase_zoom(f, f->w/2, f->h/2);
	//if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);
	if (k == '+') action_change_zoom_by_factor(f, f->w/2, f->h/2, 2);
	if (k == '-') action_change_zoom_by_factor(f, f->w/2, f->h/2, 0.5);
	if (k == 'p') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.1);
	if (k == 'm') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.1);
	if (k == 'P') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.006);
	if (k == 'M') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.006);

	if (k == 'a') action_contrast_span(f, 1/1.3);
	if (k == 'A') action_contrast_span(f, 1.3);
	//if (k == 'b') action_contrast_change(f, 1, 1);
	//if (k == 'B') action_contrast_change(f, 1, -1);
	if (k == 'n') action_qauto2(f);
	if (k == 'N') action_toggle_aqauto(f);
	if (k == 'u') action_cycle_hud(f);
	if (k == 'r') action_toggle_roi(f, x, y, m&FTR_MASK_SHIFT);
	if (k == 't') action_toggle_topography(f);
	if (k == 'c') action_toggle_p(f);
	if (k == 's') action_topography_span(f, 1/1.3);
	if (k == 'S') action_topography_span(f, 1.3);
	if (k == 'd') action_topography_Pspan(f, 1/1.3);
	if (k == 'D') action_topography_Pspan(f, 1.3);
	if (k == ',') action_screenshot(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// arrows move the viewport
	if (k > 1000 || k=='j'||k=='k'||k=='l'||k=='h') {
		int d[2] = {0, 0};
		int inc = -10;
		if (m & FTR_MASK_SHIFT  ) inc /= 10;
		if (m & FTR_MASK_CONTROL) inc *= 10;
		switch (k) {
		case 'h': case FTR_KEY_LEFT : d[0] -= inc; break;
		case 'l': case FTR_KEY_RIGHT: d[0] += inc; break;
		case 'k': case FTR_KEY_UP   : d[1] -= inc; break;
		case 'j': case FTR_KEY_DOWN : d[1] += inc; break;
		}
		if (k == FTR_KEY_PAGE_UP)   d[1] = +f->h/3;
		if (k == FTR_KEY_PAGE_DOWN) d[1] = -f->h/3;
		action_offset_viewport(f, d[0], d[1]);
	}

	if (k == '2') action_reload_image(f);
	if (k == '3' || k == ' '       ) action_flip(f, +1);
	if (k == '4' || k == FTR_KEY_BS) action_flip(f, -1);

//	// if 'k', do weird things
//	if (k == 'k') {
//		fprintf(stderr, "setting key_handler_print\n");
//		ftr_set_handler(f, "key", key_handler_print);
//	}
}


#define BAD_MIN(a,b) a<=b?a:b
#include "pickopt.c"
static char *base_name(char *p)
{
	char *b = strrchr(p, '/');
	return b ? b + 1 : p;
}
int main_cpu_single(int c, char *v[])
{
	// extract named options
	char *window_title = pick_option(&c, &v, "t", "cpu");

	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	struct pan_state e[1];
	e->i = fancy_image_open(filename_in, "r");
	e->w = e->i->w;
	e->h = e->i->h;

	// setup fonts (TODO, integrate these calls into fontu's caching stuff)
	e->font[0] = reformat_font(*xfont_4x6, UNPACKED);
	e->font[1] = reformat_font(*xfont_6x12, UNPACKED);
	e->font[2] = reformat_font(*xfont_7x13, UNPACKED);
	e->font[3] = reformat_font(*xfont_9x15, UNPACKED);
	e->font[4] = reformat_font(*xfont_10x20, UNPACKED);
	//e->font[0] = reformat_font(*xfont_5x7, UNPACKED);

	// open window
	struct FTR f = ftr_new_window(BAD_MIN(e->w,1000), BAD_MIN(e->h,800));
	ftr_change_title(&f, window_title);
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "key"   , pan_key_handler);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	//for (int i = 0; i < 5; i++) free(e->font[i].data);
	ftr_close(&f);
	//fancy_image_close(e->i);
	return r - 1;
}
int main_cpu_multi(int c, char *v[])
{
	// each input argument is an image
	// if no input arguments, read from stdin
	int n = c - 1;
	if (n > MAX_IMAGES) n = MAX_IMAGES;
	char *t[1+n];
	t[0] = "-";
	for (int i = 0; i < n; i++) t[i] = v[i+1];
	if (n == 0) n = 1;

	// read images
	struct pan_state e[1];
	e->i_num = n;
	for (int i = 0; i < n; i++) e->i_name[i] = t[i];
	for (int i = 0; i < n; i++) e->i_tab[i]  = fancy_image_open(t[i], "r");
	e->i = e->i_tab[0];
	e->w = e->i->w;
	e->h = e->i->h;

	// setup fonts (TODO, integrate these calls into fontu's caching stuff)
	e->font[0] = reformat_font(*xfont_4x6, UNPACKED);
	e->font[1] = reformat_font(*xfont_6x12, UNPACKED);
	e->font[2] = reformat_font(*xfont_7x13, UNPACKED);
	e->font[3] = reformat_font(*xfont_9x15, UNPACKED);
	e->font[4] = reformat_font(*xfont_10x20, UNPACKED);
	//e->font[0] = reformat_font(*xfont_5x7, UNPACKED);

	// open window
	struct FTR f = ftr_new_window(BAD_MIN(e->w,1000), BAD_MIN(e->h,800));
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	action_update_window_title(&f);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "key"   , pan_key_handler);
#ifdef CPU_SIGUSR
	global_f_for_sigusr = &f;
	signal(SIGUSR1, handle_signal);
	signal(SIGUSR2, handle_signal);
#endif//CPU_SIGUSR

	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	//for (int i = 0; i < 5; i++) free(e->font[i].data);
	ftr_close(&f);
	//fancy_image_close(e->i);
	return r - 1;
}

static char *help_string_name     = "cpu";
static char *help_string_version  = "cpu 2.0\n\nWritten by eml";
static char *help_string_oneliner = "interactive display of an image series";
static char *help_string_usage    = "usage:\n\tcpu image.png";
static char *help_string_long     =
"Cpu is an interface for displaying and exploring a list of images.\n"
"\n"
"This program can be used for image visualization, but this a side-effect,\n"
"just an afterthought.  The main goal of this program is to showcase\n"
"a certain philosophy for image processing: an image is an array of numbers,\n"
"it can be extrapolated towards the whole plane, and interpolated to\n"
"infinite resolution, its the values can be positive, negative, or\n"
"d-dimensional vectors, the size can be 1x1 or 10^5x10^5, and the running\n"
"time of each operation depends only on the size of the window.\n"
"\n"
"Usage: cpu image.png\n"
"   or: cat image.png | cpu\n"
"   or: cpu a.png b.png\n"
"   or: cat a.png | cpu - b.png\n"
"   or: cpu *.png\n"
"\n"
"Keys:\n"
" q,ESC  quit the program\n"
" +,-    zoom in,out by a factor of 2\n"
" a,A    increase, decrease contrast\n"
" hjkl   pan image left, down, up, right (by 10 pixels)\n"
" ^hjkl  pan image left, down, up, right (by 100 pixels)\n"
" HJKL   pan image left, down, up, right (by 1 pixel)\n"
" c      cycle between color balance modes (MIN-MAX,AVG-STD,MED-IQD)\n"
" n      normalize the contrast (globally)\n"
" N      toggle local contrast (using only the viewport data)\n"
" 2      reload the current image\n"
" SPACE  next image\n"
" BS     previous image\n"
" r      toggle display of local spectrum\n"
" u      toggle hud\n"
" t      toggle topographic mode\n"
" s,S    (in topo-mode) set vertical exaggeration\n"
" d,D    (in topo-mode) set topographic view parameter\n"
"\n"
"Mouse:\n"
" SHIFT-MOVE     center contrast under cursor\n"
" WHEEL          zoom in/out\n"
" LEFT-DRAG      pan the image domain\n"
" SHIFT-WHEEL    change local contrast span\n"
" RIGHT-CLICK    reset zoom and position\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>.\n"
;

#include "help_stuff.c"
int main_cpu(int c, char *v[])
{
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);
	return main_cpu_multi(c, v);
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_cpu(c, v); }
#endif
