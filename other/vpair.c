// gcc -std=c99 -O3 vpair.c iio.o -o vpair -lX11 -ltiff -lpng
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
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

	// 2. viewer state
	int left_w, right_w; // left_w + right_w == f->w;
	float image_scale, image_offx, image_offy;
	float hough_scale, hough_offx, hough_offy;
	int show_line;
	float line_theta, line_rho;

	// 3. pointer to the window
	struct FTR *f;
};

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	if (l < 0) l = 0;
	if (l >= pd) l = pd - 1;
	return x[pd*(j*w+i)+l];
}

static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static bool cut_line_with_segment(double out[2], double line[3],
		double p[2], double q[2])
{
	// points in "oriented" projective coordinates
	double pp[3] = {p[0], p[1], 1};
	double qq[3] = {q[0], q[1], 1};

	// sign of each point (says on which side of the line each point is)
	double sp = scalar_product(pp, line);
	double sq = scalar_product(qq, line);

	// if signs are different, the line crosses the segment
	if (sp * sq < 0) {
		// line trough points p and q
		double pq[3]; vector_product(pq, pp, qq);

		// intersection of "line" and "pq"
		double ii[3]; vector_product(ii, pq, line);

		// recover affine coordinates
		out[0] = ii[0] / ii[2];
		out[1] = ii[1] / ii[2];
		return true;
	}
	return false;
}

static bool cut_line_with_rectangle(double out_a[2], double out_b[2],
		double line[3], double rectangle[2])
{
	double w = rectangle[0];
	double h = rectangle[1];

	// four vertices of the rectangle
	double v[4][2] = { {0, 0}, {w, 0}, {w, h}, {0, h} };

	// intersections with each of the edges
	bool xP[4]; // whether it intersects
	double x[4][2]; // where it intersects
	for (int i = 0; i < 4; i++)
		xP[i] = cut_line_with_segment(x[i], line, v[i], v[ (i+1)%4 ] );

	// write output
	int n_intersections = xP[0] + xP[1] + xP[2] + xP[3];
	if (n_intersections == 2) { // generic case: 2 intersections
		int cx = 0;
		for (int i = 0; i < 4; i++)
			if (xP[i])
			{
				double *out = cx ? out_b : out_a;
				out[0] = x[i][0];
				out[1] = x[i][1];
				cx += 1;
			}
		return true;
	}
	return false;
}

static void get_line_from_polar(double line[3], double rho, double theta)
{
	line[0] = cos(theta);
	line[1] = sin(theta);
	line[2] = rho;
}

// funtion to test whether a point is inside the window
static int finsideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

// auxiliary function for drawing a red pixel
static void plot_pixel_red(int x, int y, void *e)
{
	struct FTR *f = e;
	if (finsideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}

// generic function to traverse a segment between two pixels
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (qx - px > qy - py || px - qx > qy - py) { // horizontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px + j*slope), j+py, e);
		}
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	// check consistency
	assert(f->w == e->left_w + e->right_w);

	// dump left part of the window (image)
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->left_w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->image,
				e->image_w, e->image_h, 3, i, j, l);

	// dump right part of the window (hough transform)
	for (int j = 0; j < f->h; j++)
	for (int i = e->left_w; i < f->w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->hough,
				e->hough_w, e->hough_h, 1, i-e->left_w, j, 0);

	// if necessary, show line
	if (e->show_line)
	{
		double line[3];
		get_line_from_polar(line, e->line_rho, e->line_theta);
		double rectangle[2] = {e->image_w, e->image_h};
		double p[2], q[2];
		if (cut_line_with_rectangle(p, q, line, rectangle))
			plot_segment_red(f, p[0], p[1], q[0], q[1]);
		//int i = e->line_rho;
		//int j = e->line_theta;
		//if (insideP(e->left_w, f->h, i, j))
		//{
		//	f->rgb[(j*f->w+i)*3+0] = 0;
		//	f->rgb[(j*f->w+i)*3+1] = 255;
		//	f->rgb[(j*f->w+i)*3+2] = 0;
		//}
	}
}



static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "motion F (b=%d m=%d, xy=%d %d)\n", b, m, x, y);

	struct pan_state *e = f->userdata;

	// if i'm on the left side, request line showing)
	if (x >= e->left_w && x < f->w && y >= 0 && y < f->h)
	{
		e->show_line = 1;
		e->line_theta = -y * (2*3.1416) / e->hough_h;
		e->line_rho = ((x - e->left_w) / hypot(e->hough_w, e->hough_h)) * e->hough_w;
		fprintf(stderr, "show at theta=%g rho=%g\n", e->line_theta, e->line_rho);
	} else
		e->show_line = 0;

	f->changed = 1;
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

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
	struct FTR f = ftr_new_window(e->image_w + e->hough_w, e->image_h);
	f.userdata = e;
	e->f = &f;
	e->left_w = e->image_w;
	e->right_w = e->hough_w;
	e->show_line = 0;
	f.changed = 1;

	// set handlers
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "key"   , pan_key_handler);

	// run event loop
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(e->image);
	free(e->hough);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
