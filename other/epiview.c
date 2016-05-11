// cc -O2 epiview.c iio.o -o epiview -lX11 -ltiff -lpng -lfftw3f
//
// A program for visualizing fundamental matrices
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
// * Left side: left image
// * Right side: right image
// * Both sides admit scroll (drag), zoom (wheel) and contrast changes
// * Hovering on the left shows the corresponding epipolar line on the right
// * And vice versa


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. input image pair and fundamental matrix
	int w[2], h[2], pd[2];
	float *x[2];
	double fm[9];

	// 2. view port parameters
	int octave[2];
	double zoom[2], offset[2][2];
	double aaa[2][3], bbb[2][3];

	// 3. window parameters
	int win_half;

	// 4. local state
	int scroll_domain;
	double dip_abc[3];
	int lock_transform;

	// 5. visualization details
	int head_up_display;
	struct bitmap_font font;

	// TODO: allow to click on a point to fix the visualization of its line
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

static int gmod(int x, int m)
{
	int r = x % m;
	return r < 0 ? r + m : r;
}

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (l >= pd) l = pd - 1;
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[pd*(j*w+i)+l];
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

	e->f_changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->x_octave = e->f_octave = 0;
	e->x_zoom = e->f_zoom = 1;
	e->x_offset[0] = e->x_offset[1] = e->f_offset[0] = e->f_offset[1] = 0;
	e->f_changed = 1;

	f->changed = 1;
}

static void action_change_zoom_to_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	if (F == 1) e->x_octave = 0;

	double c[2];
	window_to_image(c, e, x, y);

	e->x_zoom = 1/F;
	e->x_offset[0] = c[0] - x/e->x_zoom;
	e->x_offset[1] = c[1] - y/e->x_zoom;
	//fprintf(stderr, "\t zoom changed to %g %g {%g %g}\n", e->zoom_x, e->zoom_y, e->offset_x, e->offset_y);

	f->changed = 1;
}

static void action_change_fzoom_to_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	if (F == 1) e->f_octave = 0;

	double c[2];
	window_to_frequency(c, e, x - e->w, y);

	e->f_zoom = 1/F;
	e->f_offset[0] = c[0] - x/e->f_zoom;
	e->f_offset[1] = c[1] - y/e->f_zoom;

	f->changed = 1;
}


static void action_increase_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->x_octave < 10) {
		e->x_octave += 1;
		double fac = 1 << e->x_octave;
		if (e->x_octave < 0) fac = 1.0/(1<<-e->x_octave);
		action_change_zoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "increased octave to %d\n", e->x_octave);
}

static void action_decrease_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->x_octave > 0) {
		e->x_octave -= 1;
		double fac = 1 << e->x_octave;
		action_change_zoom_to_factor(f, x, y, fac);
	}
	else if (e->x_octave <= 0) {
		e->x_octave -= 1;
		double fac = 1.0/(1 << -e->x_octave);
		action_change_zoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "decreased octave to %d\n", e->x_octave);
}

static void action_increase_foctave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->f_octave < 10) {
		e->f_octave += 1;
		double fac = 1 << e->f_octave;
		if (e->f_octave < 0) fac = 1.0/(1<<-e->f_octave);
		action_change_fzoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "increased foctave to %d\n", e->f_octave);
}

static void action_decrease_foctave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->f_octave > 0) {
		e->f_octave -= 1;
		double fac = 1 << e->f_octave;
		action_change_fzoom_to_factor(f, x, y, fac);
	}
	else if (e->f_octave <= 0) {
		e->f_octave -= 1;
		double fac = 1.0/(1 << -e->f_octave);
		action_change_fzoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "decreased foctave to %d\n", e->f_octave);
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

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	fprintf(stderr, "center contrast at %d %d\n", x, y);
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	float c[3];
	c[0] = getsample_0(e->x, e->w, e->h, e->pd, p[0], p[1], 0);
	c[1] = getsample_0(e->x, e->w, e->h, e->pd, p[0], p[1], 1);
	c[2] = getsample_0(e->x, e->w, e->h, e->pd, p[0], p[1], 2);
	float C = (c[0] + c[1] + c[2])/3;

	e->x_b = 127.5 - e->x_a * C;

	f->changed = 1;
}

static void action_center_contrast_at_fpoint(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_frequency(p, e, x - e->x_w, y);
	float c[3];
	c[0] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 0);
	c[1] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 1);
	float C = hypot(c[0], c[1]);
	if (e->f_log)
		C = log(C);

	fprintf(stderr, "center fcontrast at %d %d (C = %g)\n", x, y, C);

	e->f_b = 127.5 - e->f_a * C;

	e->f_changed = 1;
	f->changed = 1;
}


static void action_contrast_span(struct FTR *f, float factor)
{
	fprintf(stderr, "contrast span %g\n", factor);
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->x_b)/ e->x_a;
	e->x_a *= factor;
	e->x_b = 127.5 - e->x_a * c;

	f->changed = 1;
}

static void action_contrast_fspan(struct FTR *f, float factor)
{
	fprintf(stderr, "contrast span %g\n", factor);
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->f_b)/ e->f_a;
	e->f_a *= factor;
	e->f_b = 127.5 - e->f_a * c;

	e->f_changed = 1;
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

	uint8_t fg[3] = {0, 255, 0};
	uint8_t bg[3] = {0, 0, 0};
	char buf[0x100];

	snprintf(buf, 0x100, "%g %g : %lf %lf", e->w*e->dip_a/(2*M_PI), e->h*e->dip_b/(2*M_PI), e->dip_val, e->dip_phi);
	put_string_in_rgb_image(f->rgb,f->w,f->h,e->x_w,0,fg,bg,0,&e->font,buf);
}

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	//for (int i = 0; i < f->w * f->h * 3; i++) f->rgb[i] = 0;

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
			dest[l] = float_to_byte(e->x_a * v + e->x_b);
		}
	}

	// render frequency domain
	if (e->f_changed)
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->f_w; i++)
	{
		double p[2], w[6];
		window_to_frequency(p, e, i, j);
		p[0] = gmod(p[0], e->w);
		p[1] = gmod(p[1], e->h);
		unsigned char *dest = f->rgb + 3 * (j * f->w + i + e->x_w);
		w[0] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 0);
		w[1] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 1);
		double v = hypot(w[0], w[1]);
		if (e->f_log)
			dest[0] = float_to_byte(e->f_a * log(v) + e->f_b);
		else
			dest[0] = float_to_byte(e->f_a * v + e->f_b);
		dest[1] = dest[2] = dest[0];
	}
	e->f_changed = 0;

	// render sinusoids
	if (e->show_bundle)
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->x_w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		unsigned char *dest = f->rgb + 3 * (j * f->w + i);
		double v = cos(p[0] * e->dip_a + p[1] * e->dip_b + e->dip_phi);
		//double b= getsample_0(e->x, e->w, e->h, e->pd, p[0], p[1], 0);
		//dest[0] = float_to_byte(127.5 + b*v);
		//dest[1] = dest[2] = dest[0];
		//
		dest[0] = float_to_byte(255 * (1+v)/2);
	}

	if (e->show_bundle && e->head_up_display)
		dump_hud(f);

	f->changed = 1;
}

static int symmetrize_index_inside(int i, int m)
{
	i = gmod(i, m);
	assert( i >= 0 && i < m);
	int r = 0;
	if (i >= m/2) r = i-m;
	if (i < m/2) r = i;
	return r;
}


// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "b=%d, m=%d, x=%d, y=%d\n", b, m, x, y);
	struct pan_state *e = f->userdata;

	static double ox = 0, oy = 0;

	// right side: show sinusoids
	if (m == 0 && x >= e->x_w && x < f->w && y >= 0 && y < f->h)
	{
		double p[2], v[2];
		window_to_frequency(p, e, x - e->x_w, y);
		v[0] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 0);
		v[1] = getsample_0(e->f, e->w, e->h, 2*e->pd, p[0], p[1], 1);
		double phase = atan2(v[1], v[0]);
		double norma = hypot(v[0], v[1]);
		fprintf(stderr, "f(%g %g) = %g %g : %g %g\n", p[0], p[1], v[0], v[1], phase, norma);
		e->dip_a = symmetrize_index_inside(p[0], e->w) * 2*M_PI / e->w;
		e->dip_b = symmetrize_index_inside(p[1], e->h) * 2*M_PI / e->h;
		e->dip_phi = phase;
		{ // only for HUD
			e->dip_val = norma;
			e->dip_alpha = p[0];
			e->dip_beta = p[1];
		}
		e->show_bundle = true;
		//fprintf(stderr, "show_dip %d %d (%g %g)\n", x, y, e->dip_a, e->dip_b);
		f->changed = 1;
	} else {
		e->show_bundle = false;
		f->changed = 1;
	}

	if (m == FTR_BUTTON_LEFT)
	{
		if (e->scroll_domain == 0)
			action_offset_viewport(f, x - ox, y - oy);
		if (e->scroll_domain == 1)
			action_offset_viewportf(f, x - ox, y - oy);
	}

	if (m == FTR_MASK_SHIFT  && x < e->x_w)
		action_center_contrast_at_point(f, x, y);
	if (m == FTR_MASK_SHIFT  && x >= e->x_w)
		action_center_contrast_at_fpoint(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	fprintf(stderr, "button b=%d m=%d\n", b, m);
	//if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (x < e->x_w) {
		if (b == FTR_BUTTON_UP && (m==FTR_MASK_SHIFT || m==FTR_MASK_CONTROL)) {
			action_contrast_span(f, 1/1.3); return; }
		if (b == FTR_BUTTON_DOWN && ((m==FTR_MASK_SHIFT)||m==FTR_MASK_CONTROL)){
			action_contrast_span(f, 1.3); return; }
		if (b == FTR_BUTTON_DOWN) action_increase_octave(f, x, y);
		if (b == FTR_BUTTON_UP  ) action_decrease_octave(f, x, y);
	} else {
		if (b == FTR_BUTTON_UP && (m==FTR_MASK_SHIFT || m==FTR_MASK_CONTROL)) {
			action_contrast_fspan(f, 1/1.3); return; }
		if (b == FTR_BUTTON_DOWN && ((m==FTR_MASK_SHIFT)||m==FTR_MASK_CONTROL)){
			action_contrast_fspan(f, 1.3); return; }
		if (b == FTR_BUTTON_DOWN) action_increase_foctave(f, x-e->w, y);
		if (b == FTR_BUTTON_UP  ) action_decrease_foctave(f, x-e->w, y);
	}
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
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

	if (k == '+' || k == '=') action_decrease_octave(f, f->w/2, f->h/2);
	if (k == '-') action_increase_octave(f, f->w/2, f->h/2);

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

void pan_resize(struct FTR *f, int k, int m, int x, int y)
{
	struct pan_state *e = f->userdata;
	e->xf_h = f->h;
	e->x_w = e->f_w = f->w/2;
}

#include "parsenumbers.c"
int main_pan(int c, char *v[])
{

	// process input arguments
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s matrix left.png right.png\n", *v);
		//                          0 1      2        3
		return 1;
	}
	char *matrix_str = v[1];
	char *filename_x = v[2];
	char *filename_y = v[3];

	// read images and fundamental matrix
	struct pan_state e[1];
	e->x[0] = iio_read_image_float_vec(filename_x, 0+e->w, 0+e->h, 0+e->pd);
	e->x[1] = iio_read_image_float_vec(filename_y, 1+e->w, 1+e->h, 1+e->pd);
	read_n_doubles_from_string(e->fm, matrix_str, 9);

	// init state
	e->aaa[0] = e->aaa[1] = e->aaa[2] = 1;
	e->bbb[0] = e->bbb[1] = e->bbb[2] = 1;
	e->scroll_domain = -1;
	e->show_bundle = 0;
	e->f_changed = 0;
	e->xf_h = BAD_BOUND(200, e->h, 800);
	e->x_w = e->f_w = BAD_BOUND(100, e->w, 700);
	e->contrast_mode = 0;
	e->head_up_display = true;
	e->font = *xfont9x15;
	e->font = reformat_font(e->font, UNPACKED);

	// open window
	int win_width  = BAD_BOUND(200, 2 * e->w[0], 800);
	int win_height = BAD_BOUND(200, e->h[0], 800);
	e->win_half = win_width / 2;
	struct FTR f = ftr_new_window(win_width, win_height);

	// bind state to window, initialize handlers, run loop
	f.userdata = e;
	action_reset_zoom_and_position(&f);
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

int main(int c, char *v[])
{
	return main_pan(c, v);
}
