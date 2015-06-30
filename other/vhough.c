// gcc -std=c99 -O3 vhough.c iio.o -o vhough -lX11 -ltiff -lpng
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "iio.h"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.c"

#include "xmalloc.c"

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int image_w, image_h, image_pd;
	int hough_w, hough_h, hough_pd;
	float *image, *hough;

	// 2. pointers to each window
	struct FTR *f, *g;
};

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	if (l < 0) l = 0;
	if (l >= pd) l = pd - 1;
	return x[pd*(j*w+i)+l];
}

static void pan_exposer_f(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	fprintf(stderr, "expose F\n");
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->image,
				e->image_w, e->image_h, 3, i, j, l);
}
static void pan_exposer_g(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	fprintf(stderr, "expose G\n");
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->hough,
				e->hough_w, e->hough_h, 1, i, j, 0);
}
static void pan_motion_handler_f(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "motion F (b=%d m=%d, xy=%d %d)\n", b, m, x, y);
}
static void pan_motion_handler_g(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "motion G (b=%d m=%d, xy=%d %d)\n", b, m, x, y);
}

int main_pan(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s left.png right.png\n", *v);
		//                          0 1        2
		return c;
	}
	char *filename_image  = v[1];
	char *filename_hough = v[2];

	// read images
	int w[2], h[2], pd[2];
	float *x[2];

	// read images into the "pan_state" struct
	struct pan_state e[1];
	e->image = iio_read_image_float_vec(filename_image,
			&e->image_w, &e->image_h, &e->image_pd);
	e->hough = iio_read_image_float_vec(filename_hough,
			&e->hough_w, &e->hough_h, &e->hough_pd);
	if (e->image_pd != 3 || e->hough_pd != 1)
		return fprintf(stderr, "I expect left=color & right=gray\n");

	// open windows, and cross-reference them with the state
	struct FTR f = ftr_new_window(e->image_w, e->image_h);
	struct FTR g = ftr_new_window(e->hough_w, e->hough_h);
	f.userdata = g.userdata = e;
	e->f = &f;
	e->g = &g;

	// set handlers
	ftr_set_handler(&f, "expose", pan_exposer_f);
	ftr_set_handler(&f, "motion", pan_motion_handler_f);
	ftr_set_handler(&g, "expose", pan_exposer_g);
	ftr_set_handler(&g, "motion", pan_motion_handler_g);

	// run event loop
	int r = ftr_loop_run2(&f, &g);

	// cleanup and exit (optional)
	ftr_close(&f);
	ftr_close(&g);
	free(e->image);
	free(e->hough);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
