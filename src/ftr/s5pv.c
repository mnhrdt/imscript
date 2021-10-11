// S5PV(1)                     imscript                      S5PV(1)        {{{1
//
// An explorer for Sentinel 5P Level-1 hyperspectral images
//
// cc -O3 s5pv.c iio.o ftr.o -o s5pv -lX11 -ltiff -lm -ljpeg -lpng
//
//
// Let us start small: just load an image and show three slices`
//
// TODO:
// 1. allow mouse interaction besides keyboard
// 2. add options for better local-adapted linear contrast changes
// 3. use palettes, logscale, etc
// 4. show textual info on the window, not on xterm (for easier screensharing)
// 5. load irradiance and normalize by it (toggle)
// 6. get smile-corrected slices (toggle)
// 7. project on a geographic grid (mode without slices)
// 8. flip between several bands
// 9. flip between dates
// 10. false color by combination of bands (e.g., methane-optimized ratios)

// #includes {{{1
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "ftr.h"

#include "iio.h"
#include "xmalloc.c"


// #defines {{{1

// only used for the first window at startup
#define WORSKPACE_WIDTH 1000
#define WORKSPACE_HEIGHT 1600

// useless, must be larger than 3
#define MAX_WINDOWS 10


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
	int workspace_w;
	int workspace_h;
	int nwindows;
	int window[MAX_WINDOWS][4]; // x0, y0, w, h
	float *fbuf[MAX_WINDOWS]; // temporary buffers for slices and shit
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

static int inside3P(int w, int h, int d, int x, int y, int z)
{
	return insideP(w,h,x,y) && z >= 0 && z < d;
}

// extract a (longitude, wavelength) image at the "latitude" h_index
// [this function extracts a raw sensor snapshot]
static void extract_radiance_slice_whole_0(float *s, int w, int h,
		struct svv_state *e, int h_index)
{
	assert(w == e->w);
	assert(h == e->pd);
	float (*xx)[e->w][e->pd] = (void*)e->radiance;
	float (*ss)[w] = (void*)s;
	for (int j = 0; j < e->w; j++)
	for (int i = 0; i < e->pd; i++)
		ss[i][j] = xx[h_index][j][i];
}

// extract a (longitude,latitude) image at the wavelength d_index
static void extract_radiance_slice_range_1(float *s, int w, int h,
		struct svv_state *e, int d_index, int h_offset)
{
	assert(w == e->w);
	float (*xx)[e->w][e->pd] = (void*)e->radiance;
	float (*ss)[w] = (void*)s;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int jj = e->h - (j + h_offset);
		if (inside3P(e->w,e->h,e->pd, i,jj,d_index)
				&& insideP(w,h, i,j))
			ss[j][i] = xx[jj][i][d_index];
		else
			ss[j][i] = NAN;
	}
}

// extract a (wavelength,latitude) image at the longitude w_index
static void extract_radiance_slice_range_2(float *s, int w, int h,
		struct svv_state *e, int w_index, int h_offset)
{
	assert(w == e->pd);
	float (*xx)[e->w][e->pd] = (void*)e->radiance;
	float (*ss)[w] = (void*)s;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < e->pd; i++)
	{
		int jj = e->h - (j + h_offset);
		if (inside3P(e->w,e->h,e->pd, w_index,jj,i)
				&& insideP(w,h, i,j))
			ss[j][i] = xx[jj][w_index][i];
		else
			ss[j][i] = NAN;
	}

}


// window management {{{1

static int window_add(struct svv_state *e, int x, int y, int w, int h)
{
	if (e->nwindows < MAX_WINDOWS &&
		insideP(e->workspace_w, e->workspace_h, x, y) &&
		insideP(e->workspace_w, e->workspace_h, x+w, y+h))
	{
		int r = e->nwindows;
		e->window[r][0] = x;
		e->window[r][1] = y;
		e->window[r][2] = w;
		e->window[r][3] = h;
		e->nwindows += 1;
		fprintf(stderr, "DEBUG: setting window %d: %d,%d %d,%d\n",
				r, x, y, w, h);
		return r;

	} else {
		fprintf(stderr, "WARNING: bad window (%d) %d,%d %d,%d\n",
				e->nwindows, x, y, w, h);
		return -1;
	}
}

static int window_hit(struct svv_state *e, int x, int y)
{
	for (int i = 0; i < e->nwindows; i++)
	{
		int *W = e->window[i];
		if (insideP(W[2], W[3], x - W[0], y - W[1]))
			return i;
	}
	return -1;
}

static void dump_floats_to_rgb_window(uint8_t *frgb, int fw, int fh,
		struct svv_state *e, int win, float *s, int w, int h)
{
		if (win<0 || win>=e->nwindows)
			fail("bad dump window %d\n", win);
		uint8_t (*rgb)[fw][3] = (void*)frgb;

		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			if (i >= e->window[win][2]) continue;
			if (j >= e->window[win][3]) continue;
			int ii = e->window[win][0] + i;
			int jj = e->window[win][1] + j;
			assert(insideP(fw, fh, ii, jj));

			float gg = s[j*w+i];
			if (isfinite(gg)) {
				uint8_t g = float_to_uint8(e->a * gg + e->b);
				rgb[jj][ii][0] = g;
				rgb[jj][ii][1] = g;
				rgb[jj][ii][2] = g;
			} else {
				rgb[jj][ii][0] = 100;
				rgb[jj][ii][1] = 50;
				rgb[jj][ii][2] = 0;
			}
		}
}



// state initialization {{{1

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

	// set display size
	e->workspace_w = WORSKPACE_WIDTH;
	e->workspace_h = WORKSPACE_HEIGHT;

	// reset windows
	e->nwindows = 0;
	window_add(e, 20,        20        , e->w,  e->pd              );  // 0
	window_add(e, 20,        40 + e->pd, e->w,  e->workspace_h-60-e->pd);//1
	window_add(e, 40 + e->w, 40 + e->pd, e->pd, e->workspace_h-60-e->pd);//2
	fprintf(stderr, "nwindows = %d\n", e->nwindows);
	assert(3 == e->nwindows);
	e->fbuf[0] = xmalloc(e->w  * e->pd * sizeof(float));
	e->fbuf[1] = xmalloc(e->w  * e->h * sizeof(float));
	e->fbuf[2] = xmalloc(e->pd * e->h * sizeof(float));

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

static void free_state(struct svv_state *e)
{
	for (int i = 0; i < e->nwindows; i++)
		xfree(e->fbuf[i]);
	free(e->radiance);
}

// ACTIONS {{{1

static void action_contrast_span(struct FTR *f, float factor)
{
	struct svv_state *e = f->userdata;

	float c = (127.5 - e->b)/ e->a;
	e->a *= factor;
	e->b = 127.5 - e->a * c;

	f->changed = 1;
}

static void action_offset_viewport(struct FTR *f, int dy)
//static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct svv_state *e = f->userdata;

	//e->view_offset_x -= dx/e->zoom_factor;
	//e->view_offset_y -= dy/e->zoom_factor;

	e->view_offset_y -= dy;

	f->changed = 1;
}

static void action_offset_center(struct FTR *f, int d[3])
{
	struct svv_state *e = f->userdata;

	for (int i = 0; i < 3; i++)
		e->c[i] += d[i];

	f->changed = 1;
}




// CALLBACK: svv_exposer {{{1

static void svv_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	(void)unused_x; (void)unused_y;
	//fprintf(stderr, "\n\nexpose %d %d\n", b, m);
	struct svv_state *e = f->userdata;

	// workspace: dark blue background
	for (int i = 0; i < f->w * f->h; i++)
	{
		f->rgb[3*i+0] = 0;
		f->rgb[3*i+1] = 0;
		f->rgb[3*i+2] = 100;
	}

	// window 0: background (invisible, if all goes well)
	{
		int x = e->window[0][0];
		int y = e->window[0][1];
		int w = e->window[0][2];
		int h = e->window[0][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int k = 0; k < 3; k++)
			rgb[y+j][x+i][k] = 155*(!(0==k%3)); // cyan
	}

	// window 1: background
	{
		int x = e->window[1][0];
		int y = e->window[1][1];
		int w = e->window[1][2];
		int h = e->window[1][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int k = 0; k < 3; k++)
			rgb[y+j][x+i][k] = 155*(!(1==k%3)); // magenta
	}

	// window 2: background
	{
		int x = e->window[2][0];
		int y = e->window[2][1];
		int w = e->window[2][2];
		int h = e->window[2][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int k = 0; k < 3; k++)
			rgb[y+j][x+i][k] = 155*(!(2==k%3)); // yellow
	}

	// window 0: h-slice
	{
		float *s = e->fbuf[0];
		extract_radiance_slice_whole_0(s, e->w, e->pd, e, e->c[1]);
		dump_floats_to_rgb_window(f->rgb,f->w,f->h, e,0, s,e->w,e->pd);
	}


	// window 1: d-slice
	{
		float *s = e->fbuf[1];
		int vh = e->window[1][3];
		int ho = e->view_offset_y;
		extract_radiance_slice_range_1(s, e->w, vh, e, e->c[2], ho);
		dump_floats_to_rgb_window(f->rgb,f->w,f->h, e,1, s,e->w,vh);
	}


	// window 2: w-slice
	{
		float *s = e->fbuf[2];
		int vh = e->window[2][3];
		int ho = e->view_offset_y;
		extract_radiance_slice_range_2(s, e->pd, vh, e, e->c[0], ho);
		dump_floats_to_rgb_window(f->rgb,f->w,f->h, e,2, s,e->pd,vh);
	}

	// window 1: highlight slice coords
	{
		int x = e->window[1][0];
		int y = e->window[1][1];
		int w = e->window[1][2];
		int h = e->window[1][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		int ii = e->c[0];
		if (ii >= 0 && ii < w);
			for (int j = 0; j < h; j++)
				rgb[y+j][x+ii][0] = 200;
		int jj = e->h - (e->c[1]+e->view_offset_y);
		if (jj >= 0 && jj < h)
			for (int i = 0; i < w; i++)
				rgb[y+jj][x+i][0] = 200;
	}

	// window 0: highlight slice coords
	{
		int x = e->window[0][0];
		int y = e->window[0][1];
		int w = e->window[0][2];
		int h = e->window[0][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		int ii = e->c[0];
		if (ii >= 0 && ii < w);
			for (int j = 0; j < h; j++)
				rgb[y+j][x+ii][0] = 200;
		int jj = e->c[2];
		if (jj >= 0 && jj < h)
			for (int i = 0; i < w; i++)
				rgb[y+jj][x+i][0] = 200;
	}

	// window 2: highlight slice coords
	{
		int x = e->window[2][0];
		int y = e->window[2][1];
		int w = e->window[2][2];
		int h = e->window[2][3];
		uint8_t (*rgb)[f->w][3] = (void*)f->rgb;
		int ii = e->c[2];
		if (ii >= 0 && ii < w);
			for (int j = 0; j < h; j++)
				rgb[y+j][x+ii][0] = 200;
		int jj = e->h - (e->c[1]+e->view_offset_y);
		if (jj >= 0 && jj < h)
			for (int i = 0; i < w; i++)
				rgb[y+jj][x+i][0] = 200;
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


// CALLBACK: svv_motion_handler {{{1

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

// CALLBACK: svv_key_handler {{{1
void svv_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);
	fprintf(stderr, "KEY k=%d ('%c') m=%d x=%d y=%d\n", k, k, m, x, y);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// contrast change span
	if (k == 'a') action_contrast_span(f, 1/1.3);
	if (k == 'A') action_contrast_span(f, 1.3);

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
		if (k == FTR_KEY_PAGE_UP)   d[1] = +f->h/5;
		if (k == FTR_KEY_PAGE_DOWN) d[1] = -f->h/5;
		//action_offset_viewport(f, d[0], d[1]);
		action_offset_viewport(f, d[1]);
	}

	// hjkl,np move the center of the slice trihedron
	if (k=='j'||k=='k'||k=='l'||k=='h'||k=='n'||k=='p') {
		int d[3] = {0,0,0};
		int inc = 1;
		if (k == 'j') d[1] -= inc;
		if (k == 'k') d[1] += inc;
		if (k == 'h') d[0] -= inc;
		if (k == 'l') d[0] += inc;
		if (k == 'p') d[2] -= inc;
		if (k == 'n') d[2] += inc;
		action_offset_center(f, d);
	}

}

// CALLBACK: svv_resize {{{1
static void svv_resize(struct FTR *f, int k, int m, int x, int y)
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
	struct FTR f = ftr_new_window(e->workspace_w, e->workspace_h);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , svv_key_handler);
	ftr_set_handler(&f, "button", svv_button_handler);
	//ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", svv_exposer);
	//ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	free_state(e);
	return 0;
}

int main(int c, char *v[])
{
	return main_s5pv(c, v);
}

// vim:set foldmethod=marker:
