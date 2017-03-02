// icc -std=c99 -Ofast fpantiff.c iio.o -o fpantiff -lglut -lGL -ltiff -lm
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tiff_octaves_rw.c"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.c"

#define WHEEL_FACTOR 1.4


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int w, h;
	//float *frgb;
	struct tiff_octaves t[1];

	// 2. view port parameters
	int octave;
	double zoom_factor, offset_x, offset_y;
	double a, b;

	// 3. silly options
	bool infrared;
	unsigned char *preview;
	int pw, ph;
	bool do_preview;
	int preview_position_x;
	int preview_position_y;
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

static float getsample_0(float *x, int w, int h, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[3*(j*w+i)+l];
}

static void interpolate_at(float *out, float *x, int w, int h, float p, float q)
{
	out[0] = getsample_0(x, w, h, (int)p, (int)q, 0);
	out[1] = getsample_0(x, w, h, (int)p, (int)q, 1);
	out[2] = getsample_0(x, w, h, (int)p, (int)q, 2);
}

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

// evaluate the value a position (p,q) in image coordinates
static void pixel(float *out, struct pan_state *e, double p, double q)
{
	if (p < 0 || q < 0 || p > e->t->i->w-1 || q > e->t->i->h-1) {
		int ip = p-256;
		int iq = q-256;
		int pip = mod(ip/256, 2);
		int piq = mod(iq/256, 2);
		int val = mod(pip+piq,2);
		out[0] = out[1] = out[2] = 127+val*64;
		return;
	}

	//double factor;// = 1.0 << e->octave;
	//if (e->octave >= 0)
	//	factor = 1 << e->octave;
	//else
	//	factor = 1.0 / ( 1 << - e->octave );
	

	int fmt = e->t->i->fmt;
	int bps = e->t->i->bps;
	int spp = e->t->i->spp;
	int ss = bps / 8;

	double factor = 1.0/e->zoom_factor;
	int o = e->octave;
	if (o < 0) { o = 0; factor = 1; }
	char *pix = tiff_octaves_getpixel(e->t, o, p/factor, q/factor);

	out[0] = from_sample_to_double(pix, fmt, bps);
	if (spp >= 3) {
		out[1] = from_sample_to_double(pix + ss,   fmt, bps);
		out[2] = from_sample_to_double(pix + 2*ss, fmt, bps);
		if (spp == 4 && e->infrared) {
			double o4 = from_sample_to_double(pix + 3*ss, fmt, bps);
			out[1] = 0.9*out[1] + 0.1 * o4;
		}
	} else
		out[1] = out[2] = out[0];
}

static void action_print_value_under_cursor(struct FTR *f, int x, int y)
{
	if (x<f->w && x>=0 && y<f->h && y>=0) {
		struct pan_state *e = f->userdata;
		double p[2];
		window_to_image(p, e, x, y);
		float v[3];
		pixel(v, e, p[0], p[1]);
		fprintf(stderr, "%g %g, value %g\n",p[0],p[1],v[0]);
		//float c[3];
		//interpolate_at(c, e->frgb, e->w, e->h, p[0], p[1]);
		//printf("%g\t%g\t: %g\t%g\t%g\n", p[0], p[1], c[0], c[1], c[2]);
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
	e->octave = 0;
	e->offset_x = 0;
	e->offset_y = 0;
	/*
	e->a = 1;
	e->b = 0;
	*/

	if (e->preview) {
		e->do_preview = true;
		e->zoom_factor = e->t->i->w / (double)e->pw;
		fprintf(stderr, "preview zoom factor %g\n", e->zoom_factor);
	}

	f->changed = 1;
}

static void action_exit_preview(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	e->do_preview = false;
	e->octave = 0;
	e->offset_x = x*e->zoom_factor - e->pw/2;
	e->offset_y = y*e->zoom_factor - e->ph/2;
	e->zoom_factor = 1;
	f->changed = 1;
}

static void action_contrast_change(struct FTR *f, float afac, float bshift)
{
	struct pan_state *e = f->userdata;

	e->a *= afac;
	e->b += bshift;

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

	e->a = 255 / ( M - m );
	e->b = 255 * m / ( m - M );

	f->changed = 1;
}

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	float c[3];
	pixel(c, e, p[0], p[1]);
	float C = (c[0] + c[1] + c[2])/3;

	e->b = 127.5 - e->a * C;

	f->changed = 1;
}

static void action_base_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	float c[3];
	pixel(c, e, p[0], p[1]);
	float C = (c[0] + c[1] + c[2])/3;

	e->b =  255 - e->a * C;

	f->changed = 1;
}

static void action_contrast_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->b)/ e->a;
	e->a *= factor;
	e->b = 127.5 - e->a * c;

	f->changed = 1;
}

static void action_change_zoom_to_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	if (F == 1) e->octave = 0;

	double c[2];
	window_to_image(c, e, x, y);

	e->zoom_factor = 1/F;
	e->offset_x = c[0] - x/e->zoom_factor;
	e->offset_y = c[1] - y/e->zoom_factor;
	fprintf(stderr, "\t zoom changed to %g {%g %g}\n", e->zoom_factor, e->offset_x, e->offset_y);

	f->changed = 1;
}

//static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
//{
//	struct pan_state *e = f->userdata;
//
//	double c[2];
//	window_to_image(c, e, x, y);
//
//	e->zoom_factor *= F;
//	e->offset_x = c[0] - x/e->zoom_factor;
//	e->offset_y = c[1] - y/e->zoom_factor;
//	fprintf(stderr, "\t zoom changed %g\n", e->zoom_factor);
//
//	f->changed = 1;
//}

//static void action_reset_zoom_only(struct FTR *f, int x, int y)
//{
//	struct pan_state *e = f->userdata;
//
//	action_change_zoom_by_factor(f, x, y, 1/e->zoom_factor);
//}



//static void action_increase_zoom(struct FTR *f, int x, int y)
//{
//	action_change_zoom_by_factor(f, x, y, WHEEL_FACTOR);
//}
//
//static void action_decrease_zoom(struct FTR *f, int x, int y)
//{
//	action_change_zoom_by_factor(f, x, y, 1.0/WHEEL_FACTOR);
//}


static void action_increase_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->octave < e->t->noctaves - 1) {
		e->octave += 1;
		double fac = 1 << e->octave;
		if (e->octave < 0) fac = 1.0/(1<<-e->octave);
		action_change_zoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "increased octave to %d\n", e->octave);
}

static void action_decrease_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->octave > 0) {
		e->octave -= 1;
		double fac = 1 << e->octave;
		action_change_zoom_to_factor(f, x, y, fac);
	}
	else if (e->octave <= 0) {
		e->octave -= 1;
		double fac = 1.0/(1 << -e->octave);
		action_change_zoom_to_factor(f, x, y, fac);
	}

	fprintf(stderr, "decreased octave to %d\n", e->octave);
}

static void action_toggle_infrared(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->infrared = !e->infrared;

	f->changed = 1;
}

static void dump_preview(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	if (!e->preview) return;

	f->w = e->pw;
	f->h = e->ph;
	memcpy(f->rgb, e->preview, 3*e->pw*e->ph);

	if (e->preview_position_x >= 0)
	{
		int icon_hw = (e->pw/e->zoom_factor)/2;
		int icon_hh = (e->ph/e->zoom_factor)/2;
		for (int jj = -icon_hh; jj <= icon_hh; jj++)
		for (int ii = -icon_hw; ii <= icon_hw; ii++)
		{
			int i = ii + e->preview_position_x;
			int j = jj + e->preview_position_y;
			if (i >= 0 && j >= 0 && i < f->w && j < f->h)
			{
				int idx = j * f->w + i;
				f->rgb[3*idx+2] = 255;
			}
		}
	}

	f->changed = 1;
}

static unsigned char float_to_byte(float x)
{
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

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	if (e->do_preview) {dump_preview(f); return;}

	// for every pixel in the window
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		// compute the position of this pixel in the image
		double p[2];
		window_to_image(p, e, i, j);

		// evaluate the color value of the image at this position
		float c[3];
		pixel(c, e, p[0], p[1]);

		// transform the value into RGB using the contrast change (a,b)
		unsigned char *dest = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < 3; l++)
			dest[l] = float_to_byte(e->a * c[l] + e->b);
	}
	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;
	if (e->do_preview)
	{
		e->preview_position_x = x;
		e->preview_position_y = y;
		return;
	}

	static double ox = 0, oy = 0;

	if (m == FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	if (m == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (m == FTR_MASK_SHIFT)    action_center_contrast_at_point(f, x, y);
	if (m == FTR_MASK_CONTROL)  action_base_contrast_at_point(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;
	if (e->do_preview && b == FTR_BUTTON_LEFT) {
		action_exit_preview(f, x, y); return; }

	//fprintf(stderr, "button b=%d m=%d\n", b, m);
	if (b == FTR_BUTTON_UP && (m==FTR_MASK_SHIFT || m==FTR_MASK_CONTROL)) {
		action_contrast_span(f, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && ((m==FTR_MASK_SHIFT)||m==FTR_MASK_CONTROL)){
		action_contrast_span(f, 1.3); return; }
	//if (b == FTR_BUTTON_RIGHT && m == FTR_MASK_CONTROL) {
	//	action_reset_zoom_only(f, x, y); return; }
	if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (b == FTR_BUTTON_DOWN)   action_increase_octave(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_octave(f, x, y);
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
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

	if (k == 'i') action_toggle_infrared(f);
	if (k == '+') action_decrease_octave(f, f->w/2, f->h/2);
	if (k == '-') action_increase_octave(f, f->w/2, f->h/2);
	if (k == '1') action_change_zoom_to_factor(f, x, y, 1);

	//if (k == 'p') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.1);
	//if (k == 'm') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.1);
	//if (k == 'P') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.006);
	//if (k == 'M') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.006);

	//if (k == 'a') action_contrast_change(f, 1.3, 0);
	//if (k == 'A') action_contrast_change(f, 1/1.3, 0);
	//if (k == 'b') action_contrast_change(f, 1, 1);
	//if (k == 'B') action_contrast_change(f, 1, -1);
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

// image file input/output (wrapper around iio) {{{1
#include <stdint.h>
#include "iio.h"
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

static void add_preview(struct pan_state *e, char *filename)
{
	e->preview = read_image_uint8_rgb(filename, &e->pw, &e->ph);
	e->w = e->pw;
	e->h = e->ph;
	e->do_preview = true;
}


#define BAD_MIN(a,b) a<b?a:b

#include "pickopt.c"
int main_pan(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// process input arguments
	char *filename_preview = pick_option(&c, &v, "p", "");
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s pyrpattern\n", *v);
		//                          0 1
		return 1;
	}
	char *pyrpattern = v[1];

	// read image
	struct pan_state e[1];
	int megabytes = 100;
	tiff_octaves_init(e->t, pyrpattern, megabytes);
	e->w = 700;
	e->h = 500;
	e->infrared = 4 == e->t->i->spp;
	e->preview = NULL;
	e->do_preview = false;
	e->a = 1;
	e->b = 0;
	if (*filename_preview) add_preview(e, filename_preview);

	// open window
	struct FTR f;
	if (e->preview)
		f = ftr_new_window(e->pw, e->ph);
	else
		f = ftr_new_window(BAD_MIN(e->w,700), BAD_MIN(e->h,500));
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
