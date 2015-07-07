#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

// A point in the original image domain determines a sinusoid of candidate
// positions on the Hough domain.  This function computes this sinusoid
// "rho=rho(theta)" in discretized coordinates.
static int candidate(int w, int h, int ntheta, int nrho,
		int i, int j, int i_theta,
		bool folding)
{
	float x = i - w/2;
	float y = j - h/2;
	float alpha = M_PI + atan2(y, x);
	float theta = (2 * M_PI * i_theta) / ntheta;
	if (folding) theta /= 2;
	float rho = hypot(x, y) * cos(alpha - theta);
	int i_rho = nrho   * (0.5 + rho   / hypot(w, h) );
	return i_rho;
}

// A point on the Hough domain determines a straight line in the original image.
// This function coimputes the parameters "aX+bY+c=0" of that line.
static void voterline(double abc[3], int w, int h, int ntheta, int nrho,
		int i_theta, int i_rho, int folding)
{
	float rho = hypot(w, h) * (i_rho /(float)nrho - 0.5);
	float theta = (2 * M_PI * i_theta) / ntheta;
	if (folding) theta /= 2;
	abc[0] = cos(theta);
	abc[1] = sin(theta);
	abc[2] = rho - (abc[0] * w + abc[1] * h)/2;
	fprintf(stderr, "voterline [%d %d] [%d %d] (%d %d): t=%g r=%g"
			" {%g %g %g}\n",
		w, h, ntheta, nrho, i_theta, i_rho, theta*180/M_PI, rho,
		abc[0], abc[1], abc[2]);
}

void ihough(float *transform, int ntheta, int nrho, float *imag, int w, int h,
		float minmag, bool folding)
{
	// fill counts to zero
	for (int i = 0; i < ntheta * nrho; i++)
		transform[i] = 0;

	// for each bubble
	// accumulate the corresponding Hough curve by its magnitude
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float magnitude = imag[j*w+i];

		if (magnitude < minmag)
			continue;

		for (int i_theta = 0; i_theta < ntheta; i_theta++)
		{
			int i_rho = candidate(w, h, ntheta, nrho, i, j, i_theta,
					folding);
			if (insideP(ntheta, nrho, i_theta, i_rho))
			{
				int it = i_rho * ntheta + i_theta;
				transform[it] += magnitude;
			}
		}
	}
}

#define MAIN_IHOUGH2

#ifdef MAIN_IHOUGH2
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

static void plot_transparent_red_dot(uint8_t *o, int w, int h, int x, int y)
{
	for (int j = -20; j <= 20; j++)
	for (int i = -20; i <= 20; i++)
	{
		if (hypot(i,j)>20) continue;
		int ii = x + i;
		int jj = y + j;
		float ifac = 0.99;//exp(-hypot(i,j)/10);
		if (insideP(w, h, ii, jj))
		{
			int idx = jj * w + ii;
			o[idx*4 + 0] = 128*ifac;
			o[idx*4 + 1] = 0;
			o[idx*4 + 2] = 0;
			o[idx*4 + 3] = 100;
		}
	}
	if (insideP(w, h, x, y))
	{
		int idx = y * w + x;
		o[idx*4 + 0] = 255;
		o[idx*4 + 1] = 0;
		o[idx*4 + 2] = 0;
		o[idx*4 + 3] = 255;
	}
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

// cut a line with a segment (returns true if they cut)
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

struct plot_state {
	uint8_t *x;
	int w, h;
	uint8_t rgba[4];
};

static void plot_generic_pixel(int x, int y, void *p)
{
	struct plot_state *e = p;
	if (insideP(e->w, e->h, x, y))
		for (int l = 0; l < 4; l++)
			e->x[(y*e->w+x)*4+l] = e->rgba[l];
}

static void plot_transparent_red_line(uint8_t *o, int w, int h, double l[3])
{
	double rec_from[2] = {0, 0}, rec_to[2] = {w, h}, a[2], b[2];
	bool cut = cut_line_with_rectangle(a, b, l, rec_from, rec_to);
	fprintf(stderr, "cut l[%g %g %g]:\n", l[0], l[1], l[2]);
	if (!cut) return;
	fprintf(stderr, "\t(%g %g)-(%g %g)\n", a[0], a[1], b[0], b[1]);
	struct plot_state p_red[1] = {{o, w, h, {255, 0, 0, 255}}};
	traverse_segment(a[0], a[1], b[0], b[1], plot_generic_pixel, p_red);
}


static void dump_left(char *fname_left, char *fname_right,
		int w, int h, int ntheta, int nrho,
		int click_x, int click_y, bool folding)
{
	uint8_t *alpha_left  = malloc(4 * w * h);
	uint8_t *alpha_right = malloc(4 * ntheta * nrho);

	// fill backgrounds
	for (int i = 0; i < w * h * 4; i++)
		alpha_left[i] = 0;
	for (int i = 0; i < ntheta * nrho * 4; i++)
		alpha_right[i] = 0;

	// plot transparent dot in left image
	plot_transparent_red_dot(alpha_left, w, h, click_x, click_y);

	// plot red sinusoid in right image
	for (int i_theta = 0; i_theta < ntheta; i_theta++)
	{
		int i_rho = candidate(w, h, ntheta, nrho,
				click_x, click_y, i_theta,
				folding);
		if (insideP(ntheta, nrho, i_theta, i_rho))
		{
			int idx = i_rho * ntheta + i_theta;
			alpha_right[idx*4 + 0] = 255;
			alpha_right[idx*4 + 1] = 0;
			alpha_right[idx*4 + 2] = 0;
			alpha_right[idx*4 + 3] = 255;
		}
	}

	iio_save_image_uint8_vec(fname_left, alpha_left, w, h, 4);
	iio_save_image_uint8_vec(fname_right, alpha_right, ntheta, nrho, 4);
	free(alpha_left);
	free(alpha_right);
}

static void dump_right(char *fname_left, char *fname_right,
		int w, int h, int ntheta, int nrho,
		int click_a, int click_b, bool folding)
{
	uint8_t *alpha_left  = malloc(4 * w * h);
	uint8_t *alpha_right = malloc(4 * ntheta * nrho);

	// fill backgrounds
	for (int i = 0; i < w * h * 4; i++)
		alpha_left[i] = 0;
	for (int i = 0; i < ntheta * nrho * 4; i++)
		alpha_right[i] = 0;

	// plot transparent dot in right image
	plot_transparent_red_dot(alpha_right, ntheta, nrho, click_a, click_b);

	// plot straight line in left image
	double l[3];
	voterline(l, w, h, ntheta, nrho, click_a, click_b, folding);
	plot_transparent_red_line(alpha_left, w, h, l);

	iio_save_image_uint8_vec(fname_left, alpha_left, w, h, 4);
	iio_save_image_uint8_vec(fname_right, alpha_right, ntheta, nrho, 4);
	free(alpha_left);
	free(alpha_right);
}

#include "pickopt.c"
int main(int c, char *v[])
{
	// extract optional parameters
	bool fold             =      pick_option(&c, &v, "f", NULL);
	bool omit_computation =      pick_option(&c, &v, "n", NULL);
	float click_x         = atof(pick_option(&c, &v, "x", "nan"));
	float click_y         = atof(pick_option(&c, &v, "y", "nan"));
	float click_a         = atof(pick_option(&c, &v, "a", "nan"));
	float click_b         = atof(pick_option(&c, &v, "b", "nan"));
	char *f_left          =      pick_option(&c, &v, "l", "/dev/null");
	char *f_right         =      pick_option(&c, &v, "r", "/dev/null");
	char *f_in            =      pick_option(&c, &v, "i", "-");
	char *f_out           =      pick_option(&c, &v, "o", "-");

	// process input arguments
	if (c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s ntheta nrho mtres <intensity >hough\n", *v);
		//               0  1    2      3
		return 1;
	}
	int ntheta  = atoi(v[1]);
	int nrho    = atoi(v[2]);
	float mtres = atof(v[3]);

	// read input intensities
	int w, h, pd;
	float *bubbles = iio_read_image_float_vec(f_in, &w, &h, &pd);
	if (pd != 1) return fprintf(stderr, "I expect an intensity!\n");

	if (ntheta < 0) ntheta = w;
	if (nrho < 0) nrho = h;

	// compute transform
	if (omit_computation)
		goto after_computation;
	float *transform = malloc(ntheta * nrho * sizeof*transform);
	if (!omit_computation)
		ihough(transform, ntheta, nrho, bubbles, w, h, mtres, fold);

	// save output image
	iio_save_image_float_vec(f_out, transform, ntheta, nrho, 1);

after_computation:

	// if requested, print masks
	if (isfinite(click_x + click_y))
		dump_left(f_left,f_right,w,h,ntheta,nrho,click_x,click_y,fold);
	else if (isfinite(click_a + click_b))
		dump_right(f_left,f_right,w,h,ntheta,nrho,click_a,click_b,fold);

	// cleanup and exit
	//free(bubbles);
	//free(transform);
	return 0;
}
#endif//MAIN_IHOUGH2
