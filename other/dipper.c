// gcc -std=c99 -O3 vdip.c iio.o -o vdip -lX11 -ltiff -lpng
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int image_w, image_h, image_pd;
	int hough_w, hough_h, hough_pd;
	float *image, *hough;
	double aradius;

	// 2. viewer state (images)
	int left_w, right_w; // left_w + right_w == f->w;
	float image_scale, image_offx, image_offy;
	float hough_scale, hough_offx, hough_offy;

	// 3. viewer state (decoration)
	int show_dip_bundle;
	double dip_a, dip_b;
	int dip_offset, dip_stride;
	int show_tangent;
	double tangent_space[3];
	double tangent_trans[3];

	// 4. pointer to the window
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
		double line[3], double rec_from[2], double rec_to[4])
{
	// four vertices of the rectangle
	double v[4][2] = {
		{ rec_from[0], rec_from[1] },
		{ rec_to[0]  , rec_from[1] },
		{ rec_to[0]  , rec_to[1]   },
		{ rec_from[0], rec_to[1]   }
	};

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
			if (xP[i]) {
				double *out = cx ? out_b : out_a;
				out[0] = x[i][0];
				out[1] = x[i][1];
				cx += 1;
			}
		return true;
	}
	return false;
}

static void get_line_from_polar(double line[3], double rho, double theta,
		int w, int h)
{
	line[0] = cos(theta);
	line[1] = sin(theta);
	line[2] = rho - w*cos(theta)/2 - h*sin(theta)/2;
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

static void pan_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	struct pan_state *e = f->userdata;

	// check consistency
	assert(f->w == e->left_w + e->right_w);

	// dump left part of the window (image)
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < e->left_w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->image,
				e->image_w, e->image_h, 1, i, j, l);

	// dump right part of the window (hough transform)
	for (int j = 0; j < f->h; j++)
	for (int i = e->left_w; i < f->w; i++)
	for (int l = 0; l < 3; l++)
		f->rgb[(j*f->w+i)*3+l] = getsample_0(e->hough,
				e->hough_w, e->hough_h, 1, i-e->left_w, j, 0);

	// if necessary, show dip bundle
	if (e->show_dip_bundle)
	{
		for (int i_theta = 0; i_theta < e->image_w; i_theta++)
		{
			int n_theta = e->image_w;
			int n_z = e->image_h;
			float theta = i_theta * 2 * M_PI / n_theta;
			float z = e->dip_a * cos(theta) + e->dip_b * sin(theta);
			float magic_factor = n_theta / M_PI;
			int i_z = magic_factor * z;
			for (int kk = -n_z/2 - 100;
					kk <= n_z + 100; kk += e->dip_stride)
			{
				int k = kk + i_z + e->dip_offset;
				if (insideP(f->w, f->h, i_theta, k) &&
					insideP(f->w, f->h, i_theta, k+1)&&
					insideP(f->w, f->h, i_theta, k-1))
				{
					int fidx = (k+1) * f->w + i_theta;
					f->rgb[3*fidx+0] = 255;
					f->rgb[3*fidx+1] /= 3;
					f->rgb[3*fidx+2] /= 2;
					fidx = (k+0) * f->w + i_theta;
					f->rgb[3*fidx+0] = 255;
					f->rgb[3*fidx+1] /= 3;
					f->rgb[3*fidx+2] /= 2;
					fidx = (k-1) * f->w + i_theta;
					f->rgb[3*fidx+0] = 255;
					f->rgb[3*fidx+1] /= 3;
					f->rgb[3*fidx+2] /= 2;
				}
			}
		}
	}

	// if necessary, show transformed tangent
	if (e->show_tangent)
	{
		fprintf(stderr, "spatial tangent %g %g %g\n",
				e->tangent_space[0],
				e->tangent_space[1],
				e->tangent_space[2]);
		double rec[4] = {0, 0, e->image_w, e->image_h};
		double p[2], q[2];
		if (cut_line_with_rectangle(p, q, e->tangent_space, rec, rec+2))
		{
			fprintf(stderr, "p = %g %g\n", p[0], p[1]);
			fprintf(stderr, "q = %g %g\n", q[0], q[1]);
			plot_segment_red(f, p[0], p[1], q[0], q[1]);
		}

		fprintf(stderr, "transformed tangent %g %g %g\n",
				e->tangent_trans[0],
				e->tangent_trans[1],
				e->tangent_trans[2]);
		double arad = e->aradius;
		double rek[4] = {-arad, -arad, arad, arad};
		//double p[2], q[2];
		if (cut_line_with_rectangle(p, q, e->tangent_trans, rek, rek+2))
		{
			int tside = e->hough_w;
			double alf = (tside - 1) / (2.0 * arad);
			double bet = (tside - 1) / 2.0;
			int ipa = alf * p[0] + bet + e->left_w;
			int ipb = alf * p[1] + bet;
			int iqa = alf * q[0] + bet + e->left_w;
			int iqb = alf * q[1] + bet;
			fprintf(stderr, "p = %g %g, q= %g %g\n", p[0], p[1], q[0], q[1]);
			fprintf(stderr, "ip = %d %d, iq= %d %d\n", ipa, ipb, iqa, iqb);
			plot_segment_red(f, ipa, ipb, iqa, iqb);
		}
	}
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	// right side controls
	if (x >= e->left_w && x < f->w && y >= 0 && y < f->h)
	{
		if (b == FTR_BUTTON_UP)
		{
			if (m & FTR_MASK_SHIFT)
				e->dip_stride += 1;
			else
				e->dip_offset += 1;
		}
		if (b == FTR_BUTTON_DOWN)
		{
			if (m & FTR_MASK_SHIFT)
			{
				if (e->dip_stride > 2)
					e->dip_stride -= 1;
			}
			else
				e->dip_offset -= 1;
		}
	}

	f->changed = 1;
}

static void fw_gradient_at(float g[2], float *x, int w, int h, int i, int j)
{
	double x00 = getsample_0( x, w, h, 1, i + 0, j + 0, 0);
	double x10 = getsample_0( x, w, h, 1, i + 1, j + 0, 0);
	double x01 = getsample_0( x, w, h, 1, i + 0, j + 1, 0);
	double x11 = getsample_0( x, w, h, 1, i + 1, j + 1, 0);
	g[0] = 0.5 * (x10 - x00 + x11 - x01);
	g[1] = 0.5 * (x01 - x00 + x11 - x10);
	//g[0] = x10 - x00;
	//g[1] = x01 - x00;
}


static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "motion F (b=%d m=%d, xy=%d %d)\n", b, m, x, y);

	struct pan_state *e = f->userdata;

	// if i'm on the right side, request dip bundle showing
	if (x >= e->left_w && x < f->w && y >= 0 && y < f->h)
	{
		assert(e->hough_w == e->hough_h);
		int tside = e->hough_w;
		double arad = e->aradius;
		int ia = x - e->left_w;
		int ib = y;
		e->dip_a = arad * (ia / (tside - 1.0) - 0.5);
		e->dip_b = arad * (ib / (tside - 1.0) - 0.5);
		fprintf(stderr, "vdip = %g %g\n", e->dip_a, e->dip_b);
		e->show_dip_bundle = 1;
	} else
		e->show_dip_bundle = 0;

	// if i'm on the left side, request tangent drawing
	if (x >= 0 && x < e->left_w && y >= 0 && y < f->h)
	{
		// (x,y) is the point
		// g[0], g[1] is the gradient
		float g[2];
		fw_gradient_at(g, e->image, e->image_w, e->image_h, x, y);
		e->tangent_space[0] = g[0];
		e->tangent_space[1] = g[1];
		e->tangent_space[2] = -(g[0] * x + g[1] * y);

		float theta = x * 2 * M_PI / e->image_w;
		e->tangent_trans[0] = -g[1] * sin(theta);
		e->tangent_trans[1] = g[1] * cos(theta);
		e->tangent_trans[2] = g[0];
		e->show_tangent = 1;
	} else
		e->show_tangent = 0;

	f->changed = 1;
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	//fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
	//		k, isalpha(k)?k:' ', m, x, y);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

}

static
void gradient(float *grad, float *colors, int w, int h, int pd)
{
	float *x = xmalloc(w * h * sizeof*x);

	// associated gray image (temporal to this function)
	for (int i = 0; i < w * h; i++)
		x[i] = 0;
	for (int i = 0; i < w * h; i++)
	for (int l = 0; l < pd; l++)
		x[i] += colors[i*pd+l];

	// forward-difference gradient
	for (int i = 0; i < w * h * 2; i++)
		grad[i] = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; j++)
	{
		grad[2*(j*w+i)+0] = x[j*w+i+1] - x[j*w+i];
		grad[2*(j*w+i)+1] = x[j*w+i+w] - x[j*w+i];
	}

	free(x);
}

int main_pan(int c, char *v[])
{
	// process input arguments
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s dip.png tdip.png radius\n", *v);
		//                          0 1        2       3
		return c;
	}
	char *filename_image  = v[1];
	char *filename_hough = v[2];
	double aradius = atof(v[3]);

	// read images
	int w[2], h[2], pd[2];
	float *x[2];

	// read images into the "pan_state" struct
	struct pan_state e[1];
	e->image = iio_read_image_float_vec(filename_image,
			&e->image_w, &e->image_h, &e->image_pd);
	if (0 == strcmp(filename_hough, "NONE")) {
		fail("caca");
		//e->hough_w = e->image_w;
		//e->hough_h = e->image_h;
		//e->hough = xmalloc(1 * e->hough_w * e->hough_h * sizeof(float));
		//float *g = xmalloc(2 * e->image_w * e->image_h * sizeof(float));
		//gradient(g, e->image, e->image_w, e->image_h, e->image_pd);
		//ghough(e->hough, e->hough_h
	} else
		e->hough = iio_read_image_float_vec(filename_hough,
				&e->hough_w, &e->hough_h, &e->hough_pd);
	if (e->image_pd != 1 || e->hough_pd != 1)
		return fprintf(stderr, "I expect left=gray & right=gray\n");

	// open windows, and cross-reference them with the state
	struct FTR f = ftr_new_window(e->image_w + e->hough_w, e->image_h);
	f.userdata = e;
	e->f = &f;

	// setup "state" variables
	e->left_w = e->image_w;
	e->right_w = e->hough_w;
	e->show_dip_bundle = 0;
	e->aradius = aradius;;
	e->dip_offset = 0;
	e->dip_stride = 30;
	f.changed = 1;

	// set event handlers
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
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
