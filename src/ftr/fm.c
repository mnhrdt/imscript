// c99 -Ofast ftr_mains.c ftr_x11.c iio.o -o fm -ltiff -lpng -lm -lX11 -lrt
// icc -std=c99 -Ofast ftr_mains.c ftr_x11.c iio.o -o fm -ltiff -lpng -lm -lX11 -lrt

// includes and defines {{{1
#include "seconds.c"
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "ftr.h"

#include "iio.h"

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b


// image file input/output (wrapper around iio) {{{1
static unsigned char *read_image_uint8_rgb(char *fname, int *w, int *h)
{
	int pd;
	unsigned char *x = iio_read_image_uint8_vec(fname, w, h, &pd);
	//fprintf(stderr, "riur got %d %d %d\n", *w, *h, pd);
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

static void write_image_uint8_rgb(char *fname, unsigned char *x, int w, int h)
{
	iio_write_image_uint8_vec(fname, x, w, h, 3);
}


// main_viewimage {{{1
// simple image viewer
int main_viewimage(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	free(x);
	ftr_close(&f);
	return r;
}

// main_viewimage0 {{{1
// even simpler
int main_viewimage0(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);

	sleep(4);
	return 42;
}

// main_viewimage2 {{{1
// forking image viewer
int main_viewimage2(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	ftr_fork_window_with_image_uint8_rgb(x, w, h);

	// cleanup and exit (optional)
	free(x);
	return 0;
}

// main_display {{{1
// fancier image viewer
int main_display(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	ftr_set_handler(&f, "key", ftr_handler_exit_on_ESC_or_q);
	ftr_set_handler(&f, "button", ftr_handler_stop_loop);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	free(x);
	ftr_close(&f);
	return r;
}



// main_pan {{{1

// pan: viewport data structure {{{2
#define NPYR 30

// data structure for the "panning" image viewer
//
// this data goes into the "userdata" field of the FTR window structure
//
struct pan_state {
	// the image data
	int w, h;
	unsigned char *rgb;

	// the viewing port (the actual size of the viewing port
	// is given by the "w" and "h" members of the FTR structure)
	double ntiply, offset_x, offset_y;

	// image pyramid (optional)
	unsigned char *pyr_rgb[NPYR];
	int pyr_w[NPYR], pyr_h[NPYR];
};

// change of coordinates: from window "int" pixels to image "double" point
static void window_to_image(double p[2], struct pan_state *e, int i, int j)
{
	p[0] = e->offset_x + i / e->ntiply;
	p[1] = e->offset_y + j / e->ntiply;
}

// change of coordinates: from image "double" point to window "int" pixel
static void image_to_window(int i[2], struct pan_state *e, double x, double y)
{
	i[0] = floor(x * e->ntiply - e->offset_x);
	i[1] = floor(y * e->ntiply - e->offset_y);
}





// pan: bicubic interpolation {{{2

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

static unsigned char float_to_byte(float x)
{
	int ix = x;
	if (ix < 0) return 0;
	if (ix > 255) return 255;
	return ix;
}

typedef unsigned char (*getsample_operator)(unsigned char*,int,int,int,int,int,int);

static unsigned char getsample_0(unsigned char *img, int w, int h, int pd,
		int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return img[(i+j*w)*pd + l];
}

static void bicubic_interpolation(unsigned char *result,
		unsigned char *img, int w, int h, int pd, float x, float y)
{
	x -= 1;
	y -= 1;

	getsample_operator p = getsample_0;

	int ix = floor(x);
	int iy = floor(y);
	for (int l = 0; l < pd; l++) {
		float c[4][4];
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
		float r = bicubic_interpolation_cell(c, x - ix, y - iy);
		result[l] = float_to_byte(r);
	}
}


// pan: getpixels {{{2

static void getpixel_0(unsigned char *out, unsigned char *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w) { out[0]=out[1]=out[2]=0; return; }
	if (j < 0 || j >= h) { out[0]=out[1]=out[2]=0; return; }
	out[0] = x[3*(j*w+i)+0];
	out[1] = x[3*(j*w+i)+1];
	out[2] = x[3*(j*w+i)+2];
}

static void getpixel_1(unsigned char *out, unsigned char *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	out[0] = x[3*(j*w+i)+0];
	out[1] = x[3*(j*w+i)+1];
	out[2] = x[3*(j*w+i)+2];
}

static int good_modulus(int n, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(n, -p);

	int r;
	if (n >= 0)
		r = n % p;
	else {
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}

static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p)
		r -= 1;
	if (r > p)
		r = 2*p - r;
	return r;
}

static void
getpixel_2(unsigned char *out, unsigned char *x, int w, int h, int i, int j)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	out[0] = x[3*(j*w+i)+0];
	out[1] = x[3*(j*w+i)+1];
	out[2] = x[3*(j*w+i)+2];
}

static void interpolate_at(unsigned char *out,
		unsigned char *x, int w, int h, double p, double q)
{
	getpixel_0(out, x, w, h, (int)p, (int)q);
	//bicubic_interpolation(out, x, w, h, 3, p, q);
}

static void pixel(unsigned char *out, struct pan_state *e, double p, double q)
{
	if (e->ntiply > 0.9999)
		interpolate_at(out, e->rgb, e->w, e->h, p, q);
	else {
		//if(p<0||q<0||p>=e->w||q>=e->h){out[0]=out[1]=out[2]=0;return;}
		if(p<0||q<0){out[0]=out[1]=out[2]=0;return;}
		int s = -0 - log(e->ntiply) / log(2);
		if (s < 0) s = 0;
		if (s >= NPYR) s = NPYR-1;
		int sfac = 1<<(s+1);
		int w = e->pyr_w[s];
		int h = e->pyr_h[s];
		unsigned char *rgb = e->pyr_rgb[s];
		interpolate_at(out, rgb, w, h, p/sfac, q/sfac);
	}
}

// pan: actions {{{2
static void action_print_value_at_window_position(struct FTR *f, int x, int y)
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

static void action_increment_port_offset_in_window_pixels(struct FTR *f,
		int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->ntiply;
	e->offset_y -= dy/e->ntiply;

	f->changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->ntiply = 1;
	e->offset_x = 0;
	e->offset_y = 0;

	f->changed = 1;
}

static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double c[2];
	window_to_image(c, e, x, y);

	e->ntiply *= F;
	e->offset_x = c[0] - x/e->ntiply;
	e->offset_y = c[1] - y/e->ntiply;

	f->changed = 1;
}

#define WHEELFACTOR 1.3

static void action_increase_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, WHEELFACTOR);
}

static void action_decrease_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, 1.0/WHEELFACTOR);
}



// pan: expose handler {{{2
// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "pan exposer\n");
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

// pan: motion handler {{{2
// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	static double ox = 0, oy = 0;
	if (m == 0) { // no buttons or keys, just relocate mouse center
		ox = x;
		oy = y;
	}

	if (m == FTR_BUTTON_LEFT) { // when dragging, increment offset
		action_increment_port_offset_in_window_pixels(f, x-ox, y-oy);
		ox = x;
		oy = y;
	}

	// middle button: print pixel value
	if (m == FTR_BUTTON_MIDDLE)
		action_print_value_at_window_position(f, x, y);

}




// pan: button handler {{{2
// use the mouse wheel to change zoom factor
static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	if (b == FTR_BUTTON_MIDDLE) {
		action_print_value_at_window_position(f, x, y);
		return;
	}

	if (b == FTR_BUTTON_DOWN)  action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )  action_decrease_zoom(f, x, y);
	if (b == FTR_BUTTON_RIGHT) action_reset_zoom_and_position(f);
}

// pan: key handler {{{2
void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	//fprintf(stderr, "k='%c' m=%d\n", k, m);

	if (k == '+') action_increase_zoom(f, f->w/2, f->h/2);
	if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);

	// if ESC or q, exit
	if  (k == '\033' || tolower(k)=='q')
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
		action_increment_port_offset_in_window_pixels(f, d[0], d[1]);
	}

}


// pan: multiscale pyramid {{{2
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
		a[0] = getsample_0(in, iw, ih, 3, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, 3, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, 3, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, 3, 2*i+1, 2*j+1, l);
		out[3*(ow*j + i)+l] = (a[0] + a[1] + a[2] + a[3])/4;
	}
}

static void create_pyr(struct pan_state *e)
{
	for (int s = 0; s < NPYR; s++)
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

static void free_pyr(struct pan_state *e)
{
	for (int s = 0; s < NPYR; s++)
		free(e->pyr_rgb[s]);
}

// pan: main {{{2
// panning image viewer
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
	e->offset_x = e->offset_y = 0;
	e->ntiply = 1;
	create_pyr(e);

	// open window (of constant size)
	struct FTR f = ftr_new_window(BAD_MIN(e->w,800), BAD_MIN(e->h,600));
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	//ftr_set_handler(&f, "key", ftr_handler_exit_on_ESC_or_q);
	ftr_set_handler(&f, "key", pan_key_handler);
	//int r = 0; ftr_loop_fork(&f);
	int r = ftr_loop_run(&f);


	// cleanup and exit (optional)
	ftr_close(&f);
	free(e->rgb);
	free_pyr(e);
	return r;
}

// main_twoimages {{{1
// simple image viewer (does not actually work)
int main_twoimages(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s image1 image2\n", *v);
		//                          0 1      2
		return 1;
	}
	char *filename1 = v[1];
	char *filename2 = v[2];

	// read images
	int w[2], h[2];
	unsigned char *x[2];
	x[0] = read_image_uint8_rgb(filename1, w+0, h+0);
	x[1] = read_image_uint8_rgb(filename2, w+1, h+1);

	// show image in window
	struct FTR f[2];
	f[0] = ftr_new_window_with_image_uint8_rgb(x[0], w[0], h[0]);
	f[1] = ftr_new_window_with_image_uint8_rgb(x[1], w[1], h[1]);

	sleep(4);

	// cleanup and exit (optional)
	return 0;
}

// main_twoimages2 {{{1
// forking image viewer
int main_twoimages2(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s image1 image2\n", *v);
		//                          0 1      2
		return 1;
	}
	char *filename1 = v[1];
	char *filename2 = v[2];

	// read images
	int w[2], h[2];
	unsigned char *x1 = read_image_uint8_rgb(filename1, w+0, h+0);
	unsigned char *x2 = read_image_uint8_rgb(filename2, w+1, h+1);

	// show image in window
	ftr_fork_window_with_image_uint8_rgb(x1, w[0], h[0]);
	ftr_fork_window_with_image_uint8_rgb(x2, w[1], h[1]);

	// cleanup and exit (optional)
	free(x1);
	free(x2);
	return 0;
}



// main_icrop {{{1
static void do_inline_crop_rgb(unsigned char *x, int *w, int *h, int c[4])
{
	int from[2] = {BAD_MIN(c[0], c[2]), BAD_MIN(c[1], c[3])};
	int to[2]   = {BAD_MAX(c[0], c[2]), BAD_MAX(c[1], c[3])};
	int ow = to[0] - from[0]; assert(ow > 0); assert(ow <= *w);
	int oh = to[1] - from[1]; assert(oh > 0); assert(oh <= *h);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		int o_idx = j*ow + i;
		int i_idx = (j+from[1])* *w + i+from[0];
		for (int l = 0; l < 3; l++)
			x[3*o_idx+l] = x[3*i_idx+l];
	}
	*w = ow;
	*h = oh;
}

// interactive crop
int main_icrop(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);

	// read the two corners of the crop rectangle
	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param + 0, crop_param + 1);
	ftr_wait_for_mouse_click(&f, crop_param + 2, crop_param + 3);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	write_image_uint8_rgb(filename_out, x, w, h);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

// main_icrop2 {{{1
struct icrop2_state {
	unsigned char *original_image;
	int step, ox, oy;
};

static int inbetween(int a, int b, int x)
{
	return (a <= x && x <= b) || (b <= x && x <= a);
}

void icrop2_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct icrop2_state *e = f->userdata;
	unsigned char *o = e->original_image;

	if (e->step == 0) {
		e->ox = x;
		e->oy = y;
	}

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < 3; l++)
	{
		int idx = 3*(j * f->w + i) + l;
		int inside = (e->step == 0) ?  (i>=x && j >= y) :
			inbetween(e->ox, x, i) && inbetween(e->oy, y, j);
		f->rgb[idx] = inside ? o[idx] : o[idx]/2;
	}

	f->changed = 1;
}

// interactive crop, just fancier
int main_icrop2(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image in window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	struct icrop2_state e = {x, 0, 0, 0};
	f.userdata = &e;

	// set handlers
	ftr_set_handler(&f, "motion", icrop2_motion);

	// read the two corners of the crop rectangle
	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param + 0, crop_param + 1);
	e.step = 1;
	ftr_wait_for_mouse_click(&f, crop_param + 2, crop_param + 3);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	write_image_uint8_rgb(filename_out, x, w, h);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

// main_pclick {{{1
// print to stdout the coordinates of the first mouse click
int main_pclick(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in] > pos.txt\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image on window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);

	// get the first mouse click
	int pos[2];
	ftr_wait_for_mouse_click(&f, pos+0, pos+1);

	// close the window
	ftr_close(&f);

	// print the coordinates to stdout
	printf("%d %d\n", pos[0], pos[1]);

	// cleanup and exit (optional)
	free(x);
	return 0;
}

// main_iclicks {{{1
// print to stdout an image with the given mouse clicks
int main_iclicks(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [[in] out]\n", *v);
		//                          0  1    2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	// read input image
	int w, h;
	unsigned char *x = read_image_uint8_rgb(filename_in, &w, &h);

	// show image on window
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);

	// get the first mouse click
	int pos[2], button = 0;
	while (button != FTR_BUTTON_RIGHT) {
		if (button)
			printf("%d %d\n", pos[0], pos[1]);
		ftr_wait_for_mouse_click3(&f, pos+0, pos+1, &button);
	}

	// close the window
	ftr_close(&f);

	// cleanup and exit (optional)
	free(x);
	return 0;
}


// main_random {{{1
static void draw_random(struct FTR *f, int x, int y, int k, int m)
{
	// benchmarking control
	static int cx = 0;
	cx += 1;
	if (0 == cx % 100) {
		static double os = 0;
		double s = seconds();
		double dif = s - os;
		double fps = 100/dif;
		fprintf(stdout, "CX = %d\t%g FPS\n", cx++, fps);
		os = s;
	}

	// actual drawing
	f->rgb[0] = (512341*cx)%256;
	for (int i = 1; i < f->w * f->h * 3; i++)
		f->rgb[i] = (223+f->rgb[i-1]*912359+1234567*i+54321*cx)%256;
	f->changed = 1;
}

// display an animation
int main_random(int c, char *v[])
{
	int w = 800;
	int h = 600;
	unsigned char *x = malloc(3*w*h);

	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	ftr_set_handler(&f, "idle", draw_random);
	ftr_loop_run(&f);

	ftr_close(&f);
	free(x);
	return 0;
}


// main_fire {{{1

// piecewise affine sigmoid
static float lstep(float a, float b, float t, float x)
{
	if (x < a) return 0;
	if (x > b) return t;
	return t*(x-a)/(b-a);
}


static void draw_fire(struct FTR *f, int x, int y, int k, int m)
{
	// measure time
	static int cx = 0;
	cx += 1;
	if (0 == cx % 100) {
		static double os = 0;
		double s = seconds();
		double dif = s - os;
		double fps = 100/dif;
		fprintf(stdout, "CX = %d\t%g FPS\n", cx++, fps);
		os = s;
	}

	// build palette
	static unsigned char *pal = NULL;
	if (!pal) {
		pal = malloc(3*256);
		for (int i = 0; i < 256; i++) {
			pal[3*i+0] = lstep(0,105,255,i);
			pal[3*i+1] = lstep(60,120,255,i);
			pal[3*i+2] = lstep(150,160,255,i);
		}
	}

	int num_lines_bottom = 5;
	int num_lines_hidden = 25;

	// build buffer
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
	if (!f || w != f->w || h != f->h + num_lines_hidden) { 
		w = f->w;
		h = f->h + num_lines_hidden;
		if (t) free(t);
		t = malloc(w * h * sizeof*f); 
		for (int i = 0; i < w*h; i++)
			t[i] = 104;
		cx = 0;
	}

	// draw random values at the bottom
	int p = 0;
	int rfac = cx < 75 ? 200 : 10;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = fmod(t[p] + rfac*(rand()/(1.0+RAND_MAX)),256);
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-1; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p+2*w+1] = (1.5*t[p-3*w] + 1.7 * t[p-2*w+1] 
				+ 1.5 * t[p-4*w] + 1.9 * t[p-3*w-1]
				+ 1.0 * t[p-1*w-2]
				+1.9 * t[p-4*w+1]
			) / 9.51;
	}

	// render with palette
	for (int j = 0; j < h-num_lines_hidden; j++)
	for (int i = 0; i < w; i++)
	{
		int iidx = w*(h-j-1) + i;
		int idx = (unsigned char)(lstep(105,145,255,t[iidx]));
		int pos = w*j + i;
		f->rgb[3*pos+0] = pal[3*idx+0];
		f->rgb[3*pos+1] = pal[3*idx+1];
		f->rgb[3*pos+2] = pal[3*idx+2];
	}

	f->changed = 1;
}

static void draw_minifire(struct FTR *f, int x, int y, int k, int m)
{
	// build palette
	unsigned char pal[3*256];
	for (int i = 0; i < 256; i++) {
		pal[3*i+0] = -7*i;
		pal[3*i+1] = 5*i;;
		pal[3*i+2] = 3*i;
	}

	int num_lines_bottom = 5;

	// build buffer
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
	if (!f || w != f->w || h != f->h) {
		w = f->w;
		h = f->h;
		if (t) free(t);
		t = malloc(w * h * sizeof*t);
		for (int i = 0; i < w*h; i++)
			t[i] = 104;
	}

	// draw random values at the bottom
	int p = 0;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = fmod(t[p] + 15*(rand()/(1.0+RAND_MAX)),256);
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-3; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p+2*w+1] = (1.5*t[p-3*w] + 1.7 * t[p-2*w+1]
				+ 1.5 * t[p-4*w] + 1.9 * t[p-3*w-1]
				+ 1.0 * t[p-1*w-2]
				+1.9 * t[p-4*w+1]
			) / 9.51;
	}

	// render with palette
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = (unsigned char)(t[w*(h-j-1) + i]);
		f->rgb[3*(w*j+i)+0] = pal[3*idx+0];
		f->rgb[3*(w*j+i)+1] = pal[3*idx+1];
		f->rgb[3*(w*j+i)+2] = pal[3*idx+2];
	}

	f->changed = 1;
}

static void toggle_idle2(struct FTR *f, int b, int m, int x, int y)
{
	static ftr_event_handler_t prev = NULL;
	ftr_event_handler_t new = ftr_get_handler(f, "idle");
	ftr_set_handler(f, "idle", prev);
	prev = new;
}

static void fire_resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

// display another animation
int main_fire(int c, char *v[])
{
	struct FTR f = ftr_new_window(800, 600);
	ftr_set_handler(&f, "idle", draw_fire);
	//ftr_set_handler(&f, "button", toggle_idle2);
	ftr_set_handler(&f, "button", ftr_handler_toggle_idle);
	ftr_set_handler(&f, "resize", fire_resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

// display another animation
int main_minifire(int c, char *v[])
{
	int w = 800;
	int h = 600;
	unsigned char *x = malloc(3*w*h);

	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	ftr_set_handler(&f, "idle", draw_minifire);
	//ftr_set_handler(&f, "button", ftr_handler_toggle_idle);
	ftr_set_handler(&f, "resize", fire_resize);
	ftr_loop_run(&f);

	ftr_close(&f);
	free(x);
	return 0;
}


// main_tele {{{1

static void draw_tele(struct FTR *f, int x, int y, int k, int m)
{
	// build buffer
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
	if (!f || w != f->w || h != f->h) {
		w = f->w;
		h = f->h;
		if (t) free(t);
		t = malloc(w * h * sizeof*t);
		for (int i = 0; i < w*h; i++)
			t[i] = 104;
	}

	// draw random values
	int p = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		t[p] = 255*(rand()/(1.0+RAND_MAX));
		p++;
	}

	// render
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = (unsigned char)t[w*j+i];
		f->rgb[3*(w*j+i)+0] = idx;
		f->rgb[3*(w*j+i)+1] = idx;
		f->rgb[3*(w*j+i)+2] = idx;
	}

	f->changed = 1;
}

int main_tele(int c, char *v[])
{
	int w = 320;
	int h = 200;
	unsigned char *x = malloc(3*w*h);

	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	ftr_set_handler(&f, "idle", draw_tele);
	//ftr_set_handler(&f, "button", ftr_handler_toggle_idle);
	ftr_set_handler(&f, "resize", fire_resize);
	ftr_loop_run(&f);

	ftr_close(&f);
	free(x);
	return 0;
}


// main_mandelbrot {{{1

// state of the view: crop and palette
struct mandelbrot_state {
	complex long double from, to;
	unsigned char palette[256][3];
};

// wether a number has "diverged"
static bool small(complex long double z)
{
	long double a = creall(z);
	long double b = cimagl(z);
	return a*a + b*b < 16;
}

static float mandelpoint(complex long double c, int niter)
{
	int k = 0;
	complex long double z = 0, oz = z;
	//while (cabs(z) < 4 && k++ < niter)
	//while (fabs(creall(z)) + fabs(cimagl(z)) < 4 && k++ < niter)
	while (small(z) && k++ < niter)
	{
		oz = z;
		z = z * z + c;
	}
	return cabsl(oz);//k%2;
}

struct mandel_site {
	int i, j, iter;
	complex long double c, z;
};

struct mandel_state {
	int w, h;
	complex long double from, to;
	int nsites;
	struct mandel_site *t;
	int itoff;
};

// assumes the table of sites is already allocated
static void mandel_state_start(struct mandel_state *e, int w, int h,
		complex long double from, complex long double to)
{
	long double sx = (creall(to) - creall(from)) / (w - 1);
	long double sy = (cimagl(to) - cimagl(from)) / (h - 1);
	long double ox = creall(from);
	long double oy = cimagl(from);
	fprintf(stderr, "mandel start sx = %Lg\n", sx);
	fprintf(stderr, "mandel start sy = %Lg\n", sy);
	fprintf(stderr, "mandel start ox = %Lg\n", ox);
	fprintf(stderr, "mandel start oy = %Lg\n", oy);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		struct mandel_site *s = e->t + j*w + i;
		s->c = ox + i * sx + I * (oy + j * sy);
		s->z = 0;
		s->i = i;
		s->j = j;
		s->iter = 0;
	}
	e->w = w;
	e->h = h;
	e->from = from;
	e->to = to;
	e->nsites = w * h;
	e->itoff = 0;
}

static void mandelbrot_resize(struct FTR *f, int d1, int d2, int w, int h)
{
	struct mandel_state *e = f->userdata;
	free(e->t);
	e->t = malloc(w*h*sizeof*e->t);
	mandel_state_start(e, w, h, e->from, e->to);
	memset(f->rgb, 0, 3*w*h);
}

static void action_mandel_translation(struct FTR *f, int dx, int dy)
{
	struct mandel_state *e = f->userdata;
	long double offx = dx * (creall(e->to) - creall(e->from)) / (f->w - 1);
	long double offy = dy * (cimagl(e->to) - cimagl(e->from)) / (f->h - 1);
	complex long double off = offx + I * offy;
	mandel_state_start(e, f->w, f->h, e->from + off, e->to + off);
	unsigned char *tmp = malloc(3 * f->w * f->h);
	memcpy(tmp, f->rgb, 3 * f->w * f->h);
	memset(f->rgb, 0, 3 * f->w * f->h);
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < 3; l++)
	{
		int ii = i + dx;
		int jj = j + dy; 
		if (ii < 0 || jj < 0 || ii >= f->w || jj >= f->h)
			continue;
		int idxto   = (  j * f->w +  i ) * 3 + l;
		int idxfrom = ( jj * f->w + ii ) * 3 + l;
		f->rgb[idxto] = tmp[idxfrom];
	}
	free(tmp);
}

static void mandelbrot_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == FTR_KEY_ESC || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);

	struct mandel_state *e = f->userdata;
	complex long double c = (e->from + e->to)/2;
	if (k == 'k') memset(f->rgb, 0, 3 * f->w * f->h);
	if (k == 'r') {
		mandel_state_start(e, f->w, f->h, e->from, e->to);
		memset(f->rgb, 0, 3 * f->w * f->h);
	}
	if (k == '+' || k == '=') {
		// ACTION: mandel zoom in
		mandel_state_start(e, f->w, f->h, (c+e->from)/2, (c+e->to)/2);
		unsigned char *tmp = malloc(3 * f->w * f->h);
		memcpy(tmp, f->rgb, 3 * f->w * f->h);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		for (int l = 0; l < 3; l++)
		{
			int ii = f->w/4 + i/2;
			int jj = f->h/4 + j/2;
			int idxto   = (  j * f->w +  i ) * 3 + l;
			int idxfrom = ( jj * f->w + ii ) * 3 + l;
			f->rgb[idxto] = tmp[idxfrom];
			// todo morphological erosion of black pixels
			// (to avoid block artifacts)
		}
		free(tmp);
	}
	if (k == '-') {
		// ACTION: mandel zoom out
		mandel_state_start(e, f->w, f->h, 2*e->from-c, 2*e->to-c);
		unsigned char *tmp = malloc(3 * f->w * f->h);
		memcpy(tmp, f->rgb, 3 * f->w * f->h);
		memset(f->rgb, 0, 3 * f->w * f->h);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		for (int l = 0; l < 3; l++)
		{
			int ii = 2*i - f->w/2;
			int jj = 2*j - f->h/2;
			if (ii < 0 || jj < 0 || ii >= f->w || jj >= f->h)
				continue;
			int idxto   = (  j * f->w +  i ) * 3 + l;
			int idxfrom = ( jj * f->w + ii ) * 3 + l;
			f->rgb[idxto] = tmp[idxfrom];
		}
		free(tmp);
	}
	int poff = 30;
	if (k == FTR_KEY_UP)    action_mandel_translation(f, 0, -poff);
	if (k == FTR_KEY_DOWN)  action_mandel_translation(f, 0, poff);
	if (k == FTR_KEY_LEFT)  action_mandel_translation(f, -poff, 0);
	if (k == FTR_KEY_RIGHT) action_mandel_translation(f, poff, 0);
	//if (k == FTR_KEY_UP)    action_mandel_translation(f, 0, f->h/2);
	//if (k == FTR_KEY_DOWN)  action_mandel_translation(f, 0, -f->h/2);
	//if (k == FTR_KEY_LEFT)  action_mandel_translation(f, -f->w/2, 0);
	//if (k == FTR_KEY_RIGHT) action_mandel_translation(f, f->w/2, 0);
	if (k == 'p') { e->itoff += 1; fprintf(stderr,"itoff = %d\n",e->itoff);}
	if (k == 'm') { e->itoff -= 1; fprintf(stderr,"itoff = %d\n",e->itoff);}

}

static void mandel_remove_site(struct FTR *f, struct mandel_state *e, int i)
{
	struct mandel_site *s = e->t + i;
	int k = s->j * f->w + s->i;
	//uint8_t g1 = 255.0 * cabs(s->z) / 4;
	//uint8_t g2 = 255*pow(sin(s->iter/255.0),2);
	//f->rgb[3*k+0] = g2;
	//f->rgb[3*k+1] = g2;
	//f->rgb[3*k+2] = g1;
	double idx1 = sqrt((e->itoff+s->iter)/7.0);//log(1+s->iter/30.0)-1;
	//double idx3 = log(1+s->iter/10.0);
	f->rgb[3*k+0] = 127-127*cos(idx1-3.1416/3);
	f->rgb[3*k+2] = 127-127*cos(idx1+3.1416/2);
	f->rgb[3*k+1] = 127-127*cos(idx1*1.618);
	
	e->nsites -= 1;
	struct mandel_site tmp = e->t[i];
	e->t[i] = e->t[e->nsites];
	e->t[e->nsites] = tmp;
}

static void mandel_state_one_iteration(struct FTR *f, struct mandel_state *e)
{
	for (int i = 0; i < e->nsites; i++)
	{
		struct mandel_site *s = e->t + i;
		complex long double z = s->z * s->z + s->c;
		if (!small(z))
			mandel_remove_site(f, e, i);
		else {
			s->z = z;
			s->iter += 1;
		}
	}
}

static void draw_mandelbrot(float *t, int w, int h,
		complex long double from, complex long double to, int niter)
{
	long double sx = (creall(to) - creall(from)) / (w - 1);
	long double sy = (cimagl(to) - cimagl(from)) / (h - 1);
	long double ox = creall(from);
	long double oy = cimagl(from);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		t[j*w+i] = mandelpoint(ox + i * sx + I * (oy + j * sy), niter);
}

static void draw_mandelbrot_idle(struct FTR *f, int k, int m, int x, int y)
{
	struct mandel_state *e = f->userdata;

	mandel_state_one_iteration(f, e);

	fprintf(stderr, "changed %d\n", e->nsites);
	f->changed = 1;

}

int main_mandelbrot(int c, char *v[])
{
	int w = 800;
	int h = 600;
	struct FTR f = ftr_new_window(w, h);

	struct mandel_state e[1];
	e->t = malloc(w*h*sizeof*e->t);

	complex long double from = -3 -2*I;
	complex long double to = 1+2*I;

	mandel_state_start(e, w, h, from, to);
	f.userdata = e;

	ftr_set_handler(&f, "idle", draw_mandelbrot_idle);
	ftr_set_handler(&f, "resize", mandelbrot_resize);
	ftr_set_handler(&f, "key", mandelbrot_key);
	return ftr_loop_run(&f);
}

//#define MAX_PYR 12
//
//struct msmandel_state {
//	int w, h;
//	complex long double from, to;
//
//	int pyr_levels;
//	int nsites[MAX_PYR], nactive[MAX_PYR];
//	struct mandel_site *t[MAX_PYR];
//};
//
//static void msmandel_alloc_sites(struct msmandel_state *e, int w, int h)
//{
//	e->w = w;
//	e->h = h;
//	for (int i = 0; i < MAX_PYR; i++)
//
//}



// main_events {{{1
static void print_event_key(struct FTR *f, int k, int m, int x, int y)
{
	printf("event KEY     k=%d '%c'\tm=%d (%d %d)\n", k,
			(isprint(k))?k:' ', m, x, y);
}
static void print_event_button(struct FTR *f, int k, int m, int x, int y)
{
	printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
}
static void print_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
}
static void print_event_resize(struct FTR *f, int b, int m, int x, int y)
{
	printf("event RESIZE  b=%d\tm=%d (%d %d)\n", b, m, x, y);
}

int main_events(int c, char *v[])
{
	struct FTR f = ftr_new_window(320, 200);
	for (int i = 0; i < 3 * f.w * f.h; i++)
		f.rgb[i] = 0xa0*!(i%3);
	f.changed = 1;
	fprintf(stderr, "i'm here!\n");
	ftr_set_handler(&f, "key", print_event_key);
	ftr_set_handler(&f, "button", print_event_button);
	ftr_set_handler(&f, "motion", print_event_motion);
	ftr_set_handler(&f, "resize", print_event_resize);
	fprintf(stderr, "i'm still here!\n");
	return ftr_loop_run(&f);
}

// main_paint {{{1
static void paint_event_button(struct FTR *f, int b, int m, int x, int y)
{
	if (x < f->w && y < f->h && x >= 0 && y >= 0)
	{
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = f->rgb[3*idx+1] = f->rgb[3*idx+2] = 255;
	}
	f->changed = 1;
}
static void paint_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	if (x < f->w && y < f->h && x >= 0 && y >= 0)
	{
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = f->rgb[3*idx+1] = f->rgb[3*idx+2] = 127;
	}
	f->changed = 1;
}

int main_paint(int c, char *v[])
{
	struct FTR f = ftr_new_window(320, 200);
	ftr_set_handler(&f, "button", paint_event_button);
	ftr_set_handler(&f, "motion", paint_event_motion);
	return ftr_loop_run(&f);
}

// main_mini {{{1
int main_mini(int c, char *v[])
{
	struct FTR f = ftr_new_window(320, 200);
	//f.rgb[1+3*(f.w*10+10)] = 255;
	return ftr_loop_run(&f);
}
int main_mini0(int c, char *v[])
{
	struct FTR f = ftr_new_window(320, 200);
	return 0;
}
int main_mini1(int c, char *v[])
{
	struct FTR f = ftr_new_window(320, 200);
	sleep(3);
	return 0;
}

// main {{{1

#define MAIN(x) { #x, main_ ## x }

static const struct { char *n; int(*f)(int,char*[]); } mains[] = {
	MAIN(mini),
	MAIN(mini0),
	MAIN(mini1),
	MAIN(viewimage),
	MAIN(viewimage0),
	MAIN(viewimage0),
	MAIN(display),
	MAIN(pan),
	MAIN(twoimages),
	MAIN(twoimages2),
	MAIN(icrop),
	MAIN(icrop2),
	MAIN(pclick),
	MAIN(iclicks),
	MAIN(random),
	MAIN(fire),
	MAIN(minifire),
	MAIN(tele),
	MAIN(mandelbrot),
	MAIN(events),
	MAIN(paint),
	{ 0, 0 }
};

static char *base_name(char *p)
{
	char *b = strrchr(p, '/');
	return b ? b + 1 : p;
}

int main(int c, char *v[])
{
	char *t, *s = base_name(*v);
	if (s[0] == 'f' && s[1] == 'm')
		return main(c - 1, v + 1);
	int i = 0;
	while ((t = mains[i++].n))
		if (0 == strcmp(s, t))
			return mains[i-1].f(c, v);
	return 1;
}

// vim:set foldmethod=marker:
