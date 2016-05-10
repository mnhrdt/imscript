// gcc -std=c99 -g dosdo.c iio.o -o dosdo -lX11 -ltiff -lpng -lfftw3f
//
// A program for visualizing the spatial and frequential domain of an image
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>

#include <stdint.h>
#include "iio.h"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.c"

#define WHEEL_FACTOR 1.4

#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "xfont9x15.c"

// VISUALIZATION:
// * A single window split in half.
// * Left side: viewport into the image
// * Right side: viewport into its fourier transform
// * Both sides admit scroll (drag) and zoom (wheel)
// * Hovering on the right shows the corresponding basis function on the left.


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. input image pair
	int w, h, pd;
	float *x, *f; // TODO: add DCT, DHT (DST ?)

	// 2. view port parameters
	int x_octave, f_octave;
	double x_zoom, f_zoom, x_offset[2], f_offset[2];
	double x_a, x_b, f_a, f_b;
	int f_log;

	// 3. window parameters
	int x_w, f_w, xf_h;

	// 4. local state
	int scroll_domain;

	// 5. visualization details
	int contrast_mode;
	int head_up_display;
};

// change of coordinates: from window "int" pixels to image "double" point
static void window_to_image(double p[2], struct pan_state *e, int i, int j)
{
	p[0] = e->x_offset[0] + i / e->x_zoom;
	p[1] = e->x_offset[1] + j / e->x_zoom;
}

// change of coordinates: from window "int" pixels to image "double" point
static void window_to_frequency(double p[2], struct pan_state *e, int i, int j)
{
	p[0] = e->f_offset[0] + i / e->f_zoom;
	p[1] = e->f_offset[1] + j / e->f_zoom;
}

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[pd*(j*w+i)+l];
}

//static void interpolate_at(float *out, float *x, int w, int h, float p, float q)
//{
//	out[0] = getsample_0(x, w, h, (int)p, (int)q, 0);
//	out[1] = getsample_0(x, w, h, (int)p, (int)q, 1);
//	out[2] = getsample_0(x, w, h, (int)p, (int)q, 2);
//}

static int mod(int n, int p)
{
	if (p < 1)
		fail("bad modulus %d\n", p);

	int r;
	if (n >= 0)
		r = n % p;
	else {
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	if (r < 0 || r >= p)
		fprintf(stderr, "mod(%d,%d)=%d\n",n,p,r);
	assert(r >= 0);
	assert(r < p);
	return r;
}


#include<stdarg.h>
static void img_debug(float *x, int w, int h, int pd, const char *fmt, ...)
{
	return;
	va_list ap;
	char fname[FILENAME_MAX];
	va_start(ap, fmt);
	vsnprintf(fname, FILENAME_MAX, fmt, ap);
	va_end(ap);
	fprintf(stderr, "IMG_DEBUG(%dx%d,%d) \"%s\"\n", w, h, pd, fname);
	iio_save_image_float_vec(fname, x, w, h, pd);
}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->x_offset[0] -= dx/e->x_zoom;
	e->x_offset[1] -= dy/e->x_zoom;

	f->changed = 1;
}

static void action_offset_viewportf(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->f_offset[0] -= dx/e->f_zoom;
	e->f_offset[1] -= dy/e->f_zoom;

	f->changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->x_octave = e->f_octave = 0;
	e->x_zoom = e->f_zoom = 1;
	e->x_offset[0] = e->x_offset[1] = e->f_offset[0] = e->f_offset[1] = 0;

	f->changed = 1;
}

static void action_contrast_change(struct FTR *f, float afac, float bshift)
{
	struct pan_state *e = f->userdata;

	e->x_a *= afac;
	e->x_b += bshift;

	f->changed = 1;
}

static void action_qauto(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	//float m = INFINITY, M = -m;
	float m = 0, M = 255;
	//int pid = 3;
	//for (int i = 0; i < 3 * e->pyr_w[pid] * e->pyr_h[pid]; i++)
	//{
	//	float g = e->pyr_rgb[pid][i];
	//	m = fmin(m, g);
	//	M = fmax(M, g);
	//}

	e->x_a = 255 / ( M - m );
	e->x_b = 255 * m / ( m - M );

	f->changed = 1;
}

static void action_contrast_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->x_b)/ e->x_a;
	e->x_a *= factor;
	e->x_b = 127.5 - e->x_a * c;

	f->changed = 1;
}

static void action_save_shot(struct FTR *f)
{
	static int shot_counter = 1;
	char fname[FILENAME_MAX];
	snprintf(fname, FILENAME_MAX, "/tmp/vnav_shot_%d.png", shot_counter);
	iio_save_image_uint8_vec(fname, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "saved shot \"%s\"\n", fname);
	shot_counter += 1;
}


static void action_toggle_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->head_up_display = !e->head_up_display;
	f->changed = 1;
}

static unsigned char float_to_byte(float x)

{
	if (isnan(x)) return 0;
	if (x < 0) return 0;
	if (x > 255) return 255;
	// set gamma=2
	//float r = x * x / 255;
	//
	//float n = x / 255;
	//float r = (n*n)*n;
	//return r*255;
	return x;
}


static void dump_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	// Data to put on HUD:
	float first_row, last_row, dip_a, dip_b;
	double p[2], q[2];
	window_to_image(p, e, 0, 0);
	window_to_image(q, e, 0, e->xf_h - 1);
	first_row = p[1];
	last_row = q[1];

	uint8_t fg[3] = {0, 255, 0};
	uint8_t bg[3] = {0, 0, 0};
	char buf[0x100];

	snprintf(buf, 0x100, "row: %d", (int)first_row);
	//put_string_in_rgb_image(f->rgb,f->w,f->h,0,0,fg,bg,0,&e->font, buf);
}

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	for (int i = 0; i < f->w * f->h * 3; i++) f->rgb[i] = 0;

	// render spatial domain
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->x_w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		unsigned char *dest = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < 3; l++)
		{
			float v = getsample_0(e->x, e->w, e->h, e->pd,
					p[0], p[1], l);
			dest[l] = float_to_byte(v);
		}
	}

	// render frequency domain
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->f_w; i++)
	{
		double p[2];
		window_to_frequency(p, e, i, j);
		unsigned char *dest = f->rgb + 3 * (j * f->w + i + e->x_w);
		for (int l = 0; l < 3; l++)
		{
			float v = getsample_0(e->f, e->w, e->h, 2*e->pd,
					p[0], p[1], l);
			dest[l] = float_to_byte(v);
		}
	}

	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	static double ox = 0, oy = 0;

	// right side: show sinusoids
	if (m == 0 && m >= e->x_w && x < f->w && y >= 0 && y < f->h)
	{
		//double arad = e->aradius;
		//int ia = x - e->strip_w;
		//int ib = y;
		//e->dip_a = arad * (ia / (tside - 1.0) - 0.5);
		//e->dip_b = arad * (ib / (tside - 1.0) - 0.5);
		//e->show_dip_bundle = true;
		//fprintf(stderr, "show_dip %d %d (%g %g)\n", x, y, e->dip_a, e->dip_b);
		f->changed = 1;
	} else {
		//e->show_dip_bundle = false;
		f->changed = 1;
	}

	if (m == FTR_BUTTON_LEFT)
	{
		if (e->scroll_domain == 0)
			action_offset_viewport(f, x - ox, y - oy);
		if (e->scroll_domain == 1)
			action_offset_viewportf(f, x - ox, y - oy);
	}


	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	fprintf(stderr, "button b=%d m=%d\n", b, m);
	if (b == FTR_BUTTON_UP && (m==FTR_MASK_SHIFT || m==FTR_MASK_CONTROL)) {
		action_contrast_span(f, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && ((m==FTR_MASK_SHIFT)||m==FTR_MASK_CONTROL)){
		action_contrast_span(f, 1.3); return; }
	//if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	//if (b == FTR_BUTTON_DOWN)   action_increase_octave(f, x, y);
	//if (b == FTR_BUTTON_UP  )   action_decrease_octave(f, x, y);
	//if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
	if (b == FTR_BUTTON_LEFT) e->scroll_domain = x >= e->x_w;
}

void key_handler_print(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "key pressed %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isprint(k)?k:' ', m, x, y);

	//if (k == '+' || k == '=') action_decrease_octave(f, f->w/2, f->h/2);
	//if (k == '-') action_increase_octave(f, f->w/2, f->h/2);

	if (k == 'u') action_toggle_hud(f);

	if (k == ',') action_save_shot(f);

	if (k == 'n') action_qauto(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

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
}

#define BAD_BOUND(a,x,b) ((a)<(x)?((x)<(b)?(x):(b)):(a))

int main_pan(int c, char *v[])
{

	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s img fft\n", *v);
		//                          0 1   2
		return 1;
	}
	char *filename_x = v[1];
	char *filename_f = v[2];

	// read images
	struct pan_state e[1];
	int ww, hh, ppdd;
	e->x = iio_read_image_float_vec(filename_x, &e->w, &e->h, &e->pd);
	e->f = iio_read_image_float_vec(filename_f, &ww, &hh, &ppdd);
	if (e->w != ww || e->h != hh || 2*e->pd != ppdd)
		fail("domain size mismatch (%d %d %d) != (%d %d %d)\n",
				e->w, e->h, e->pd, ww, hh, ppdd);

	// init state
	e->x_a = 1;
	e->x_b = 0;
	e->f_a = 10;
	e->f_b = 0;
	e->f_log = 1;
	e->contrast_mode = 0;
	e->head_up_display = true;
	e->scroll_domain = -1;
	e->xf_h = BAD_BOUND(200, e->h, 800);
	e->x_w = e->f_w = BAD_BOUND(200, e->w, 1000);

	// open window
	struct FTR f = ftr_new_window(e->x_w + e->f_w, e->xf_h);
	f.userdata = e;
	action_reset_zoom_and_position(&f);
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
