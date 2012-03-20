// given a photo of a picture and four given points,
// produce a rectangular crop of the image whithin the given quadrilater

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "xmalloc.c"
#include "getpixel.c"

#include "cmphomod.c"

static int quadrant_signature(float midx, float midy, float p[4][2])
{
	int q[4];
	for (int i = 0; i < 4; i++)
		q[i] = 2*(p[i][1] > midy) + (p[i][0] > midx);
	return q[0] + 2*q[1] + 4*q[2] + 8*q[3];
}


// Sort four points so that they are ordered like that:
//
// p2 | p3
// ---+---
// p0 | p1
//
static bool canonicalize_point_ordering_inplace(float x[4][2])
{
	// copy input
	float p[4][2];
	for (int i = 0; i < 8; i++)
		p[i/2][i%2] = x[i/2][i%2];

	for (int i = 0; i < 4; i++)
		fprintf(stderr, "p[%d] = %g %g\n", i, p[i][0], p[i][1]);

	// find max and min on each dimension
	int imx = 0, imy = 0, iMx = 0, iMy = 0;
	for (int i = 1; i < 4; i++) {
		if (p[i][0] < p[imx][0]) imx = i;
		if (p[i][1] < p[imy][1]) imy = i;
		if (p[i][0] > p[iMx][0]) iMx = i;
		if (p[i][1] > p[iMy][1]) iMy = i;
	}

	// find two perpendicular lines separating each of the two non-extrema
	float midx = 0, midy = 0;
	for (int i = 0; i < 4; i++) {
		if (i != imx && i != iMx)
			midx += p[i][0]/2;
		if (i != imy && i != iMy)
			midy += p[i][1]/2;
	}
	fprintf(stderr, "midx = %g\n", midx);
	fprintf(stderr, "midy = %g\n", midy);

	// find the quadrant for each point
	int quadrant[4] = {-1, -1, -1, -1};
	int qs = quadrant_signature(midx, midy, p);
	for (int i = 0; i < 4; i++) {
		int qidx = 2*(p[i][1] > midy) + (p[i][0] > midx);
		assert(qidx >= 0);
		assert(qidx < 4);
		if (quadrant[qidx] > 0)
		{
			assert(qs != 15);
			switch(qs) {
			case 9:
				break;
			case 6:
				break;
			default:
				return false;
			}
			//x[0][0] = p[imx][0]; x[0][1] = p[imy][1];
			//x[1][0] = p[iMx][0]; x[1][1] = p[imy][1];
			//x[2][0] = p[imx][0]; x[2][1] = p[iMy][1];
			//x[3][0] = p[iMx][0]; x[3][1] = p[iMy][1];
			//return false;
		}
		quadrant[qidx] = i;
		fprintf(stderr, "quadrant[%d] = %d\n", qidx, i);
	}

	// save the re-ordered points
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
			x[i][j]  = p[quadrant[i]][j];

	return true;
}

// Build a rectangle that approximates the given four points.
// The input points should be ordered like that:
//
// p2 | p3
// ---+---
// p0 | p1
//
static void compute_rectangular_fit(float np[4][2], float p[4][2])
{
	assert(p[0][0] < p[1][0]);
	assert(p[2][0] < p[3][0]);
	assert(p[0][1] < p[2][1]);
	assert(p[1][1] < p[3][1]);

	float xleft  = floor((p[0][0] + p[2][0])/2);
	float xright = ceil ((p[1][0] + p[3][0])/2);
	float ylow   = floor((p[0][1] + p[1][1])/2);
	float yhigh  = ceil ((p[2][1] + p[3][1])/2);
	assert(xleft < xright);
	assert(ylow < yhigh);
	float xspan = xright - xleft;
	float yspan = yhigh - ylow;
	fprintf(stderr, "xleft = %g\n", xleft);
	fprintf(stderr, "xright = %g\n", xright);
	fprintf(stderr, "ylow = %g\n", ylow);
	fprintf(stderr, "yhigh = %g\n", yhigh);
	fprintf(stderr, "xspan = %g\n", xspan);
	fprintf(stderr, "yspan = %g\n", yspan);

	np[0][0] = 0;     np[0][1] = 0;
	np[1][0] = xspan; np[1][1] = 0;
	np[2][0] = 0;     np[2][1] = yspan;
	np[3][0] = xspan; np[3][1] = yspan;
}


// compute the homography given by the images of four points
static void compute_homography_from_point_pairs(double H[3][3],
		float from[4][2], float to[4][2])
{
	double f[4][2], t[4][2];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++) {
			f[i][j] = from[i][j];
			t[i][j] = to[i][j];
		}

	homography_from_4corresp(f[0], f[1], f[2], f[3],
	                         t[0], t[1], t[2], t[3], H);
}

// apply an homography to a point
static void apply_homography(float y[2], float x[2], double *H)
{
	float z[3];
	z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
	z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
	z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = z[0]/z[2];
	y[1] = z[1]/z[2];
}

// fill-in a vector field determined by a global homographic transform
static void fill_homographic_flow_field(float *ff, int w, int h, double H[3][3])
{
	float (*f)[w][2] = (void*)ff;
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++) {
			float p[2] = {i, j}, q[2];
			apply_homography(q, p, H[0]);
			for (int l = 0; l < 2; l++)
				f[j][i][l] = q[l] - p[l];
		}
}


// interpolate a cell, bilinearly
static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

// interpolate a cell, by nearest neighbor
static float interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	// return a;
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

// interpolate a cell, by an arbitrary method
static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return interpolate_nearest(a, b, c, d, x, y);
	//case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	default: fail("caca de vaca");
	}
}

#include "bicubic.c"

// interpolate an image at a given sub-pixelic point
static void general_interpolate(float *result,
		float *x, int w, int h, int pd, float p, float q,
		int m) // method
{
	if (m == 3) {
		bicubic_interpolation(result, x, w, h, pd, p, q);
	} else {
		int ip = floor(p);
		int iq = floor(q);
		for (int l = 0; l < pd; l++) {
			float a = getsample_0(x, w, h, pd, ip  , iq  , l);
			float b = getsample_0(x, w, h, pd, ip  , iq+1, l);
			float c = getsample_0(x, w, h, pd, ip+1, iq  , l);
			float d = getsample_0(x, w, h, pd, ip+1, iq+1, l);
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			result[l] = v;
		}
	}
}


// pull back an image by a given vector field
static void pull_back(float *yy, int yw, int yh, float *ff,
		float *xx, int xw, int xh, int pd)
{
	float (*y)[yw][pd] = (void*)yy;
	float (*f)[yw][2] = (void*)ff;
	for (int j = 0; j < yh; j++)
		for (int i = 0; i < yw; i++) {
			float p[2] = {i, j};
			float q[2] = {i + f[j][i][0], j + f[j][i][1]};
			float val[pd];
			general_interpolate(val, xx, xw, xh, pd, q[0], q[1], 3);
			for (int l = 0; l < pd; l++)
				y[j][i][l] = val[l];
		}
}


// Crop a rectangle from an image.  The rectangle is given by 4 points, which
// need not form a rectangle, but after the crop they are deformed into a
// rectangle.
//
// Warning! this function re-orders the input points
static void deframe(float *y, int *out_w, int *out_h,
		float *x, int in_w, int in_h, int pd, float points[4][2])
{
	bool good_data = canonicalize_point_ordering_inplace(points);
	//if (!good_data)
	//	fail("the four given points are too far from a rectangle");

	float cpoints[4][2];
	compute_rectangular_fit(cpoints, points);
	assert(cpoints[0][0] == 0);
	assert(cpoints[0][1] == 0);
	assert(cpoints[3][0] == floor(cpoints[3][0]));
	assert(cpoints[3][1] == floor(cpoints[3][1]));

	double H[3][3];
	compute_homography_from_point_pairs(H, cpoints, points);

	fprintf(stderr, "invH =");
	for (int i = 0; i < 9; i++)
		fprintf(stderr, " %g", H[0][i]);
	fprintf(stderr, "\n");

	*out_w = cpoints[3][0];
	*out_h = cpoints[3][1];
	fprintf(stderr, "out_w = %d\n", *out_w);
	fprintf(stderr, "out_h = %d\n", *out_h);
	assert(*out_w < in_w);
	assert(*out_h < in_h);
	float *f = xmalloc(*out_w * *out_h * 2 * sizeof*f);

	fill_homographic_flow_field(f, *out_w, *out_h, H);

	pull_back(y, *out_w, *out_h, f, x, in_w, in_h, pd);

	free(f);

}

#include "iio.h"
int main(int c, char *v[])
{
	if (c != 10 && c != 11 && c != 12) {
		fprintf(stderr, "usage:\n\t"
				"%s ax ay bx by cx cy dx dy [in [out]]\n", *v);
	//                       0  1  2  3  4  5  6  7  8   9   10
		return EXIT_FAILURE;;
	}
	float points[4][2];
	for (int i = 0; i < 8; i++)
		points[i/2][i%2] = atof(v[1+i]);
	char *filename_in  = c > 9 ? v[9] : "-";
	char *filename_out = c > 10 ? v[10] : "-";
	int w[2], h[2], pd;
	float *x = iio_read_image_float_vec(filename_in, w, h, &pd);
	float *y = xmalloc(w[0] * h[0] * pd * sizeof*y);

	deframe(y, w+1, h+1, x, *w, *h, pd, points);

	iio_save_image_float_vec(filename_out, y, w[1], h[1], pd);
	free(y);
	free(x);
	return EXIT_SUCCESS;
}
