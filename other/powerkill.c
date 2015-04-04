// icc -std=c99 -Ofast powerkill.c iio.o -o powerkill -lX11 -ltiff -ljpeg -lpng -lz -lm
#include <assert.h>
#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "iio.h"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.c"

#include "xmalloc.c"

#define WHEEL_FACTOR 1.4
#define MAX_PYRAMID_LEVELS 30

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int w, h, pd;
	complex float *fft;
	uint8_t *mask;

	// 2. view port parameters
	double zoom_factor, offset_x, offset_y;
	double a, b, log_offset;

	// 3. image pyramid
	complex float *pyr_fft[MAX_PYRAMID_LEVELS];
	int pyr_w[MAX_PYRAMID_LEVELS], pyr_h[MAX_PYRAMID_LEVELS];

	// 4. dragging states
	bool painting;
	int paint_handle[4];
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
	i[0] = floor(x * e->zoom_factor - e->offset_x);
	i[1] = floor(y * e->zoom_factor - e->offset_y);
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
		r = p - (-n) % p;
		if ((int)r == p)
			r = 0;
	}
	return r;
}

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	if (l < 0) l = 0;
	if (l >= pd) l = pd - 1;
	return x[pd*(j*w+i)+l];
}

static uint8_t bygetpixel_s(uint8_t *xx, int w, int h, int i, int j)
{
	uint8_t (*x)[w] = (void*)xx;
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	return x[j][i];
}

static uint8_t bygetpixel_0(uint8_t *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[j*w+i];
}

static float getsample_s(float *xx, int w, int h, int pd, int i, int j, int l)
{
	float (*x)[w][pd] = (void*)xx;
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	l = good_modulus(l, pd);
	return x[j][i][l];
}

static void interpolate_at(complex float *out,
		complex float *x, int w, int h, int pd, float p, float q)
{
	float *fout = (float*)out;
	for (int i = 0; i < pd; i++)
	{
		fout[2*i+0] = getsample_s((float*)x, w, h, 2*pd, p, q, 2*i+0);
		fout[2*i+1] = getsample_s((float*)x, w, h, 2*pd, p, q, 2*i+1);
	}
}

// evaluate the value a position (p,q) in image coordinates
static void cpixel(complex float *out, struct pan_state *e, double p, double q)
{
	if (e->zoom_factor > 0.9999)
		interpolate_at(out, e->fft, e->w, e->h, e->pd, p, q);
	else {
		//static int first_run = 1;
		//if (first_run) {
		//	fprintf(stderr, "create pyramid\n");
		//	void create_pyramid(struct pan_state *e);
		//	create_pyramid(e);
		//	first_run = 0;
		//}
		//if(p<0||q<0){out[0]=out[1]=out[2]=0;return;}
		int s = -0 - log(e->zoom_factor) / log(2);
		if (s < 0) s = 0;
		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
		int sfac = 1<<(s+1);
		int w = e->pyr_w[s];
		int h = e->pyr_h[s];
		complex float *fft = e->pyr_fft[s];
		interpolate_at(out, fft, w, h, e->pd, p/sfac, q/sfac);
	}
	int ip = p;
	int iq = q;
	uint8_t m = bygetpixel_s(e->mask, e->w, e->h, ip, iq);
	if (!m)
		for (int l = 0; l < e->pd; l++)
			out[l] = 0;
}

//static void action_print_value_under_cursor(struct FTR *f, int x, int y)
//{
//	if (x<f->w && x>=0 && y<f->h && y>=0) {
//		struct pan_state *e = f->userdata;
//		double p[2];
//		window_to_image(p, e, x, y);
//		float c[3];
//		interpolate_at(c, e->frgb, e->w, e->h, p[0], p[1]);
//		printf("%g\t%g\t: %g\t%g\t%g\n", p[0], p[1], c[0], c[1], c[2]);
//	}
//}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->zoom_factor;
	e->offset_y -= dy/e->zoom_factor;

	f->changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->zoom_factor = 1;
	e->offset_x = e->w/2;
	e->offset_y = e->h/2;
	e->a = 30;
	e->b = -150;
	e->log_offset = 1;

	e->painting = false;

	f->changed = 1;
}

static void action_contrast_change(struct FTR *f, float afac, float bshift)
{
	struct pan_state *e = f->userdata;

	e->a *= afac;
	e->b += bshift;

	f->changed = 1;
}

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	complex float c[e->pd];
	cpixel(c, e, p[0], p[1]);
	//float C = (c[0] + c[1] + c[2])/3;
	double X = log(e->log_offset + cabs(*c));
	e->b = 127.5 - e->a * X;
	f->changed = 1;
}

static void action_clear_mask(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	for (int i = 0; i < e->w * e->h; i++)
		e->mask[i] = 1;
	f->changed = 1;
}

static void action_symmetrize_mask(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
		if (!e->mask[j*e->w + i])
		{
			int ii = e->w - i;
			int jj = e->h - j;
			e->mask[jj*e->w + ii] = 0;
		}
	f->changed = 1;
}

static void action_radialize_mask(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	int rtab[e->w+e->h];
	for (int i = 0; i < e->w+e->h; i++)
		rtab[i] = 1;

	// build table of "bad" radii
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
		if (!e->mask[j*e->w + i])
		{
			int ii = i; if (ii > e->w/2) ii = i - e->w;
			int jj = j; if (jj > e->h/2) jj = j - e->h;
			int r = lrint(hypot(ii, jj));
			rtab[r] = 0;
		}

	// re-fill mask, discarding bad radii
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
	{
		int ii = i; if (ii > e->w/2) ii = i - e->w;
		int jj = j; if (jj > e->h/2) jj = j - e->h;
		int r = lrint(hypot(ii, jj));
		if (!rtab[r])
			e->mask[j*e->w+i] = 0;
	}

	f->changed = 1;
}

//static void action_qauto(struct FTR *f)
//{
//	struct pan_state *e = f->userdata;
//
//	float m = INFINITY, M = -m;
//	int pid = 3;
//	for (int i = 0; i < 3 * e->pyr_w[pid] * e->pyr_h[pid]; i++)
//	{
//		float g = e->pyr_fft[pid][i];
//		m = fmin(m, g);
//		M = fmax(M, g);
//	}
//
//	e->a = 255 / ( M - m );
//	e->b = 255 * m / ( m - M );
//
//	f->changed = 1;
//}

static void action_contrast_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->b)/ e->a;
	e->a *= factor;
	e->b = 127.5 - e->a * c;

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

static void action_kill_windowed_rectangle(struct pan_state *e, int x, int y)
{
	int f_good[2] = { BAD_MIN(e->paint_handle[0], x),
			  BAD_MIN(e->paint_handle[1], y) };
	int t_good[2] = { BAD_MAX(e->paint_handle[0], x),
			  BAD_MAX(e->paint_handle[1], y) };
	double from[2], to[2];
	window_to_image(from, e, f_good[0], f_good[1]);
	window_to_image(to,   e, t_good[0], t_good[1]);
	int ifrom[2], ito[2];
	ifrom[0] = good_modulus(from[0], e->w);
	ifrom[1] = good_modulus(from[1], e->h);
	ito[0] = good_modulus(to[0], e->w);
	ito[1] = good_modulus(to[1], e->h);
	fprintf(stderr, "\twill kill %g %g to %g %g\n", from[0], from[1], to[0], to[1]);
	fprintf(stderr, "\twill kill %d %d to %d %d\n", ifrom[0], ifrom[1], ito[0], ito[1]);
	int rw = lrint(fabs(from[0] - to[0]));
	int rh = lrint(fabs(from[1] - to[1]));
	for (int j = 0; j <= rh; j++)
	for (int i = 0; i <= rw; i++)
	{
		int ii = good_modulus(from[0] + i, e->w);
		int jj = good_modulus(from[1] + j, e->h);
		e->mask[jj*e->w+ii] = 0;
	}
}

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

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		complex float c[e->pd];
		cpixel(c, e, p[0], p[1]);
		unsigned char *cc = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < e->pd; l++) {
			double X = log(e->log_offset + cabs(c[l]));
			double Y = e->a * X + e->b;
			cc[l] = float_to_uint8(Y);
		}
		if (e->pd < 3)
			cc[1] = cc[2] = cc[0];
		int ip = lrint(p[0]);
		int iq = lrint(p[1]);
		if (insideP(e->w, e->h, ip, iq) && !e->mask[e->w*iq+ip])
			cc[0] = cc[1] = cc[2] = 0;

	}
	if (e->painting) {
		int from[2] = {
			BAD_MIN(e->paint_handle[0], e->paint_handle[2]),
			BAD_MIN(e->paint_handle[1], e->paint_handle[3])
		};
		int to[2] = {
			BAD_MAX(e->paint_handle[0], e->paint_handle[2]),
			BAD_MAX(e->paint_handle[1], e->paint_handle[3])
		};
		fprintf(stderr, "expose while painting %d %d => %d %d\n",
				from[0], from[1], to[0], to[1]);
		for (int j = from[1]; j <= to[1]; j++)
		for (int i = from[0]; i <= to[0]; i++)
			if (insideP(f->w, f->h, i, j))
			{
				uint8_t *cc = f->rgb + 3 * (j * f->w + i);
				cc[0] = cc[1] = cc[2] = 0;
			}
	}
	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	static double ox = 0, oy = 0;

	if (m == FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	//if (m == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (m == FTR_MASK_SHIFT)    action_center_contrast_at_point(f, x, y);
	//
	struct pan_state *e = f->userdata;
	if (e->painting) {
		e->paint_handle[2] = x;
		e->paint_handle[3] = y;
		f->changed = 1;
	}

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "button b=%d m=%d\n", b, m);
	if (b == FTR_BUTTON_UP && m == FTR_MASK_SHIFT) {
		action_contrast_span(f, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && m == FTR_MASK_SHIFT) {
		action_contrast_span(f, 1.3); return; }

	struct pan_state *e = f->userdata;

	if (b == FTR_BUTTON_RIGHT)
	{
		e->paint_handle[0] = x;
		e->paint_handle[1] = y;
		e->painting = true;
		fprintf(stderr, "paint handle from %d %d\n", x, y);
	}

	if (e->painting && b == -FTR_BUTTON_RIGHT)
	{
		e->painting = false;
		fprintf(stderr, "end painting at %d %d\n", x, y);
		action_kill_windowed_rectangle(e, x, y);
		f->changed = 1;
	}

	//if (b == FTR_BUTTON_RIGHT && m == FTR_MASK_CONTROL) {
	//	action_reset_zoom_only(f, x, y); return; }
	//if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	//if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
}

void key_handler_print(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "key pressed %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);

	if (k == 'c') action_clear_mask(f);
	if (k == 's') action_symmetrize_mask(f);
	if (k == 'r') action_radialize_mask(f);

	//if (k == '+') action_increase_zoom(f, f->w/2, f->h/2);
	//if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);
	if (k == '+') action_change_zoom_by_factor(f, f->w/2, f->h/2, 2);
	if (k == '=') action_change_zoom_by_factor(f, f->w/2, f->h/2, 2);
	if (k == '-') action_change_zoom_by_factor(f, f->w/2, f->h/2, 0.5);
	if (k == 'p') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.1);
	if (k == 'm') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.1);
	if (k == 'P') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.006);
	if (k == 'M') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.006);

	//if (k == 'a') action_contrast_change(f, 1.3, 0);
	//if (k == 'A') action_contrast_change(f, 1/1.3, 0);
	//if (k == 'b') action_contrast_change(f, 1, 1);
	//if (k == 'B') action_contrast_change(f, 1, -1);
	//if (k == 'n') action_qauto(f);
	if (k == '.')  action_reset_zoom_and_position(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
	{
		struct pan_state *e = f->userdata;
		if (e->painting)
			e->painting = 0;
		else
			ftr_notify_the_desire_to_stop_this_loop(f, 1);
	}

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

	// if 'k', do weird things
	if (k == 'k') {
		fprintf(stderr, "setting key_handler_print\n");
		ftr_set_handler(f, "key", key_handler_print);
	}
}

static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < pd; l++)
	{
		float a[4];
		a[0] = getsample_0(in, iw, ih, pd, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, pd, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, pd, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, pd, 2*i+1, 2*j+1, l);
		out[pd*(ow*j + i)+l] = (a[0] + a[1] + a[2] + a[3])/4;
	}
}

// type of a "zoom-out" function
typedef void (*zoom_out_function_t)(float*,int,int,float*,int,int,int);

static void create_pyramid(struct pan_state *e)
{
	zoom_out_function_t z = zoom_out_by_factor_two;
	int pd = e->pd;
	for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
	{
		int      lw   = s ? e->pyr_w  [s-1] : e->w   ;
		int      lh   = s ? e->pyr_h  [s-1] : e->h   ;
		complex float   *lfft = s ? e->pyr_fft[s-1] : e->fft;
		int      sw   = ceil(lw / 2.0);
		int      sh   = ceil(lh / 2.0);
		complex float *sfft = xmalloc(sw * sh * pd * sizeof*sfft);
		z((float*)sfft, sw, sh, (float*)lfft, lw, lh, 2*pd);
		e->pyr_w[s]   = sw;
		e->pyr_h[s]   = sh;
		e->pyr_fft[s] = sfft;
	}
}

#include "pickopt.c"

int main_pan(int c, char *v[])
{
	// process input arguments
	char *imask_option = pick_option(&c, &v, "m", "");
	if (c > 4) {
		fprintf(stderr, "usage:\n\t"
				"%s [in.fft [out.fft [out.mask]]]\n", *v);
		//                0  1       2        3
		return c;
	}
	char *filename_in   = c > 1 ? v[1] : "-";
	char *filename_out  = c > 2 ? v[2] : "-";
	char *filename_mask = c > 3 ? v[3] : NULL;
	char *filename_imask = *imask_option ? imask_option : NULL;


	// read image
	struct pan_state e[1];
	int pd;
	e->fft = (void*)iio_read_image_float_vec(filename_in, &e->w, &e->h,&pd);
	if (pd != 2 && pd != 6)
		return fprintf(stderr, "input must be a fft (got pd=%d)\n", pd);
	e->pd = pd / 2;
	create_pyramid(e);
	if (filename_imask) {
		int mw, mh;
		e->mask = iio_read_image_uint8(filename_imask,&mw,&mh);
		if (mw != e->w || mh != e->h)
			return fprintf(stderr, "input mask bad size\n");
		fprintf(stderr, "using input mask from file \"%s\"\n",
				filename_imask);
	} else {
		e->mask = malloc(e->w * e->h);
		memset(e->mask, 1, e->w * e->h);
	}

	// open window
	struct FTR f = ftr_new_window(BAD_MIN(e->w,512), BAD_MIN(e->h,512));
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "key"   , pan_key_handler);
	int r = ftr_loop_run(&f);

	// apply computed mask
	for (int i = 0; i < e->w * e->h; i++)
		if (!e->mask[i])
			for (int l = 0; l < e->pd; l++)
				e->fft[i*e->pd+l] = 0;

	// save output images
	iio_save_image_float_vec(filename_out, (void*)e->fft, e->w,e->h, pd);
	if (filename_mask)
		iio_save_image_uint8_vec(filename_mask, e->mask, e->w, e->h, 1);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
