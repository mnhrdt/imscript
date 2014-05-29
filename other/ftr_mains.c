// c99 -Ofast ftr_mains.c ftr_x11.c iio.o -o fm -ltiff -lpng -lm -lX11 -lrt
// icc -std=c99 -Ofast ftr_mains.c ftr_x11.c iio.o -o fm -ltiff -lpng -lm -lX11 -lrt

#include "seconds.c"
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "ftr.h"

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

static void write_image_uint8_rgb(char *fname, unsigned char *x, int w, int h)
{
	iio_save_image_uint8_vec(fname, x, w, h, 3);
}


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



#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b

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




typedef unsigned char (*getsample_operator)(unsigned char*,int,int,int,int,int,int);

static unsigned char getsample_0(unsigned char *img, int w, int h, int pd,
		int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return img[(i+j*w)*pd + l];
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

static unsigned char float_to_byte(float x)
{
	int ix = x;
	if (ix < 0) return 0;
	if (ix > 255) return 255;
	return ix;
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
	//assert(r >= 0);
	//assert(r < p);
	return r;
}

static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p)
		r -= 1;
	if (r > p)
		r = 2*p - r;
	//assert(r >= 0);
	//assert(r < p);
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
	//bicubic_interpolation(out, x, w, h, 3, p, q);
	getpixel_0(out, x, w, h, (int)p, (int)q);
}

static inline void pixel(unsigned char *out, struct pan_state *e, double p, double q)
{
	if (e->ntiply >= 1)
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

#define FTR_BUTTON_LEFT 1
#define FTR_BUTTON_MIDDLE 2
#define FTR_BUTTON_RIGHT 3
#define FTR_BUTTON_UP 4
#define FTR_BUTTON_DOWN 5

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
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	static double ox = 0, oy = 0;
	if (m == 0) { // just relocate mouse center
		ox = x;
		oy = y;
	}
	if (m == 256) { // when dragging, increment offset
		action_increment_port_offset_in_window_pixels(f, x-ox, y-oy);
		ox = x;
		oy = y;
		f->changed = 1;
	}

	// middle button: print pixel value
	//if (b == FTR_BUTTON_MIDDLE) {
	if (m == 512) 
		action_print_value_at_window_position(f, x, y);

}


static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double c[2];
	window_to_image(c, e, x, y);

	e->ntiply *= F;
	e->offset_x = c[0] - x/e->ntiply;
	e->offset_y = c[1] - y/e->ntiply;
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

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->ntiply = 1;
	e->offset_x = 0;
	e->offset_y = 0;
}




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
	f->changed = 1;
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
	ftr_set_handler(&f, "key", ftr_handler_exit_on_ESC_or_q);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(e->rgb);
	free_pyr(e);
	return r;
}

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



static void do_inline_crop_rgb(unsigned char *x, int *w, int *h, int crop[4])
{
	int from[2] = {BAD_MIN(crop[0], crop[2]), BAD_MIN(crop[1], crop[3])};
	int to[2]   = {BAD_MAX(crop[0], crop[2]), BAD_MAX(crop[1], crop[3])};
	int ow = to[0] - from[0]; //assert(ow > 0); assert(ow <= *w);
	int oh = to[1] - from[1]; //assert(oh > 0); assert(oh <= *h);
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
	ftr_wait_for_mouse_click(&f, crop_param+0, crop_param+1, NULL, NULL);
	ftr_wait_for_mouse_click(&f, crop_param+2, crop_param+3, NULL, NULL);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	write_image_uint8_rgb(filename_out, x, w, h);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

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
	ftr_wait_for_mouse_click(&f, crop_param+0, crop_param+1, NULL, NULL);
	e.step = 1;
	ftr_wait_for_mouse_click(&f, crop_param+2, crop_param+3, NULL, NULL);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	write_image_uint8_rgb(filename_out, x, w, h);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

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
	ftr_wait_for_mouse_click(&f, pos+0, pos+1, NULL, NULL);

	// close the window
	ftr_close(&f);

	// print the coordinates to stdout
	printf("%d %d\n", pos[0], pos[1]);

	// cleanup and exit (optional)
	free(x);
	return 0;
}


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
		f->changed = 1;
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
	ftr_set_handler(&f, "key", print_event_key);
	ftr_set_handler(&f, "button", print_event_button);
	ftr_set_handler(&f, "motion", print_event_motion);
	ftr_set_handler(&f, "resize", print_event_resize);
	return ftr_loop_run(&f);
}

int main(int c, char *v[])
{
	int (*f)(int,char*[]);
	if (c < 2) return fprintf(stderr, "name a main\n");
	else if (0 == strcmp(v[1], "viewimage"))  f = main_viewimage;
	else if (0 == strcmp(v[1], "viewimage0")) f = main_viewimage0;
	else if (0 == strcmp(v[1], "viewimage2")) f = main_viewimage2;
	else if (0 == strcmp(v[1], "display"))    f = main_display;
	else if (0 == strcmp(v[1], "pan"))        f = main_pan;
	else if (0 == strcmp(v[1], "twoimages"))  f = main_twoimages;
	else if (0 == strcmp(v[1], "twoimages2")) f = main_twoimages2;
	else if (0 == strcmp(v[1], "icrop"))      f = main_icrop;
	else if (0 == strcmp(v[1], "icrop2"))     f = main_icrop2;
	else if (0 == strcmp(v[1], "pclick"))     f = main_pclick;
	else if (0 == strcmp(v[1], "random"))     f = main_random;
	else if (0 == strcmp(v[1], "fire"))       f = main_fire;
	else if (0 == strcmp(v[1], "minifire"))   f = main_minifire;
	else if (0 == strcmp(v[1], "events"))     f = main_events;
	//else if (0 == strcmp(*v, "simplest"))     f = main_simplest;
	//else if (0 == strcmp(*v, "simplest2"))    f = main_simplest2;
	else return fprintf(stderr, "bad main \"%s\"\n", v[1]);
	return f(c-1, v+1);
}
