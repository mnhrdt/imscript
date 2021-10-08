// S5PV(1)                     imscript                      S5PV(1)        {{{1
//
// An explorer for Sentinel 5P Level-1 hyperspectral images
//
// cc -O3 s5pv.c iio.o ftr.o -o s5pv -lX11 -ltiff -lm -ljpeg -lpng
//
//
// Let us start small: just load an image and show three slices`

// #includes {{{1
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h> // for getpid only


#include "ftr.h"

#include "iio.h"
#include "xmalloc.c"


// #defines {{{1

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define EARTH_RADIUS 6378000.0
#define WHEEL_FACTOR 2.0


// svv_state {{{1

// this goes inside the "userdata" field of the 
struct svv_state {
	// 0. image data (so far, radiance-only)
	int w, h, pd;
	float *radiance;

	// 1. view port parameters
	double a, b; // linear contrast change
	double view_offset_y;

	// 2. interaction parameters
	int c[3]; // center of the slice trihedron

	// 3. display shit, subwindows, etc
	// ...
};


// utility functions {{{1

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



// state initialization functions {{{1

static void init_state(struct svv_state *e, char *filename_in)
{
	// load image data
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	fprintf(stderr, "image of size = %dx%d,%d\n", w, h, pd);

	// set image data inside the state (just to begin with)
	e->w = w;
	e->h = h;
	e->pd = pd;
	e->radiance = x;

	// set viewport
	e->a = 1e8;
	e->b = 0;
	e->view_offset_y = 300;

	// set slice trihedron center
	e->c[0] = w/2;
	e->c[1] = h/2;
	e->c[2] = pd/2;

	// not used: save the reat data plus central slices
	//
	//iio_write_image_float_vec("/tmp/full_radiance.npy", x, w, h, pd);
	//// save the three central slides
	//float *y = malloc((w+h+pd)*(w+h+pd)*sizeof*y);
	//float (*xx)[w][pd] = (void*)x;

	//{
	//	float (*yy)[w] = (void*)y;
	//	for (int j = 0; j < h; j++)
	//	for (int i = 0; i < w; i++)
	//		yy[j][i] = xx[j][i][pd/2];
	//	iio_write_image_float("/tmp/slyce_pd.npy", y, w, h);
	//}

	//{
	//	float (*yy)[pd] = (void*)y;
	//	for (int j = 0; j < h; j++)
	//	for (int i = 0; i < pd; i++)
	//		yy[j][i] = xx[j][w/2][i];
	//	iio_write_image_float("/tmp/slyce_w.npy", y, pd, h);
	//}

	//{
	//	float (*yy)[pd] = (void*)y;
	//	for (int j = 0; j < w; j++)
	//	for (int i = 0; i < pd; i++)
	//		yy[j][i] = xx[h/2][j][i];
	//	iio_write_image_float("/tmp/slyce_h.npy", y, pd, w);
	//}

	//free(y);
}







// CALLBACK: svv_exposer {{{1

static void svv_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	(void)unused_x; (void)unused_y;
	//fprintf(stderr, "\n\nexpose %d %d\n", b, m);
	struct svv_state *e = f->userdata;

	// dark blue background
	for (int i = 0; i < f->w * f->h; i++)
	{
		f->rgb[3*i+0] = 0;
		f->rgb[3*i+1] = 0;
		f->rgb[3*i+2] = 100;
	}

	// just one slice shit
	float (*xx)[e->w][e->pd] = (void*)e->radiance;
	for (int j = 0; j < e->w; j++)
	for (int i = 0; i < e->pd; i++)
	{
		float gg = xx[e->h/2][j][i];
		uint8_t g = float_to_uint8(e->a * gg + e->b);
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		rgb[20+j][20+i][0] = g;
		rgb[20+j][20+i][1] = g;
		rgb[20+j][20+i][2] = g;
	}
}


// CALLBACK: svv_button_handler {{{1
static void svv_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "BUT b=%d m=%d x=%d y=%d\n", b, m, x, y);
	//if (b == FTR_BUTTON_UP && m & FTR_MASK_CONTROL) {
	//	action_cycle_view(f, +1, x, y); return; }
	//if (b == FTR_BUTTON_DOWN && m & FTR_MASK_CONTROL) {
	//	action_cycle_view(f, -1, x, y); return; }
	//if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	//if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	//if (b == FTR_BUTTON_MIDDLE) action_print_data_under_cursor(f, x, y);
}


// CALLBACK: pan_motion_handler {{{1

// update offset variables by dragging
static void svv_motion_handler(struct FTR *f, int unused_b, int m, int x, int y)
{
	//(void)unused_b;
	//struct pan_state *e = f->userdata;
	//static double ox = 0, oy = 0;
	//if (m & FTR_BUTTON_LEFT) action_offset_viewport(f, x - ox, y - oy);
	//ox = x;
	//oy = y;

	//if (e->show_vertdir)
	//{
	//	e->vdx = x;
	//	e->vdy = y;
	//	f->changed = 1;
	//}
}

// CALLBACK: pan_key_handler {{{1
void svv_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

}

// CALLBACK: pan_resize {{{1
static void pan_resize(struct FTR *f, int k, int m, int x, int y)
{
	//request_repaints(f);
}


// main {{{1
int main_s5pv(int c, char *v[])
{
	// process input arguments
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s l1_XXXXXX.nc,radiance\n", *v);
		//                          0 1
		return c;
	}
	char *filename_in = v[1];


	// start state
	struct svv_state e[1];
	init_state(e, filename_in);



	// open window
	struct FTR f = ftr_new_window(1000, 1000);
	f.userdata = e;
	f.changed = 1;
	//ftr_set_handler(&f, "key"   , pan_key_handler);
	//ftr_set_handler(&f, "button", pan_button_handler);
	//ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", svv_exposer);
	//ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return 0;
}

int main(int c, char *v[])
{
	return main_s5pv(c, v);
}

// vim:set foldmethod=marker:
