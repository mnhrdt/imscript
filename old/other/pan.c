#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "ftr_mini.h"
#include "iio.h"

#define WHEEL_FACTOR 1.3
#define MAX_PYRAMID_LEVELS 30

// image file input/output (wrapper around iio) {{{1
static unsigned char *read_image_uint8_rgb(char *fname, int *w, int *h)
{
	int pd;
	unsigned char *x = iio_read_image_uint8_vec(fname, w, h, &pd);
	if (pd == 3) return x;
	unsigned char *y = malloc(3**w**h);
	for (int i = 0; i < *w**h; i++) {
		switch(pd) {
		case 1:
			y[3*i+0] = y[3*i+1] = y[3*i+2] = x[i];
			break;
		case 2:
			y[3*i+0] = x[2*i+0];
			y[3*i+1] = y[3*i+2] = x[2*i+1];
			break;
		default:
			assert(pd > 3);
			for (int l = 0; l < 3; l++)
				y[3*i+l] = x[pd*i+l];
			break;
		}
	}
	free(x);
	return y;
}

// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int w, h;
	unsigned char *rgb;

	// 2. view port parameters
	double zoom_factor, offset_x, offset_y;

	// 3. image pyramid
	unsigned char *pyr_rgb[MAX_PYRAMID_LEVELS];
	int pyr_w[MAX_PYRAMID_LEVELS], pyr_h[MAX_PYRAMID_LEVELS];
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

static unsigned char getsample_0(unsigned char *x, int w, int h,
						int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[3*(j*w+i)+l];
}

static void interpolate_at(unsigned char *out,
		unsigned char *x, int w, int h, double p, double q)
{
	out[0] = getsample_0(x, w, h, (int)p, (int)q, 0);
	out[1] = getsample_0(x, w, h, (int)p, (int)q, 1);
	out[2] = getsample_0(x, w, h, (int)p, (int)q, 2);
}

// evaluate the value a position (p,q) in image coordinates
static void pixel(unsigned char *out, struct pan_state *e, double p, double q)
{
	if (e->zoom_factor > 0.9999)
		interpolate_at(out, e->rgb, e->w, e->h, p, q);
	else {
		if(p<0||q<0){out[0]=out[1]=out[2]=0;return;}
		int s = -0 - log(e->zoom_factor) / log(2);
		if (s < 0) s = 0;
		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
		int sfac = 1<<(s+1);
		int w = e->pyr_w[s];
		int h = e->pyr_h[s];
		unsigned char *rgb = e->pyr_rgb[s];
		interpolate_at(out, rgb, w, h, p/sfac, q/sfac);
	}
}

static void action_print_value_under_cursor(struct FTR *f, int x, int y)
{
	if (x<f->w && x>=0 && y<f->h && y>=0) {
		struct pan_state *e = f->userdata;
		double p[2];
		window_to_image(p, e, x, y);
		unsigned char c[3];
		interpolate_at(c, e->rgb, e->w, e->h, p[0], p[1]);
		printf("%g\t%g\t: %d\t%d\t%d\n", p[0], p[1], c[0], c[1], c[2]);
	}
}

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
	e->offset_x = 0;
	e->offset_y = 0;

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

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		unsigned char *c = f->rgb + 3 * (j * f->w + i);
		pixel(c, e, p[0], p[1]);
	}
	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	static double ox = 0, oy = 0;

	if (m == FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	if (m == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (k == '+') action_increase_zoom(f, f->w/2, f->h/2);
	if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);

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
		action_offset_viewport(f, d[0], d[1]);
	}
}

static void zoom_out_by_factor_two(unsigned char *out, int ow, int oh,
		unsigned char *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < 3; l++)
	{
		float a[4];
		a[0] = getsample_0(in, iw, ih, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, 2*i+1, 2*j+1, l);
		out[3*(ow*j + i)+l] = (a[0] + a[1] + a[2] + a[3])/4;
	}
}

static void create_pyramid(struct pan_state *e)
{
	for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
	{
		int            lw   = s ? e->pyr_w  [s-1] : e->w  ;
		int            lh   = s ? e->pyr_h  [s-1] : e->h  ;
		unsigned char *lrgb = s ? e->pyr_rgb[s-1] : e->rgb;
		int            sw   = ceil(lw / 2.0);
		int            sh   = ceil(lh / 2.0);
		unsigned char *srgb = malloc(3 * sw * sh);
		zoom_out_by_factor_two(srgb, sw, sh, lrgb, lw, lh);
		e->pyr_w[s]   = sw;
		e->pyr_h[s]   = sh;
		e->pyr_rgb[s] = srgb;
	}
}

static void free_pyramid(struct pan_state *e)
{
	for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
		free(e->pyr_rgb[s]);
}

#define BAD_MIN(a,b) a<b?a:b

int main_pan(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	struct pan_state e[1];
	e->rgb = read_image_uint8_rgb(filename_in, &e->w, &e->h);
	create_pyramid(e);

	// open window
	struct FTR f = ftr_new_window(BAD_MIN(e->w,800), BAD_MIN(e->h,600));
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "key", pan_key_handler);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(e->rgb);
	free_pyramid(e);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
