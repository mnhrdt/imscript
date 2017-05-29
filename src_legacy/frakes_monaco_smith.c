// frakes-monaco-smith turbulence code
//
// 1. increase image resolution fourfold using bilinear interpolation
// 2. inverse Laplacian emphasis (PARAMETER sigma)
// 3. CGI
// 	3.1. PARAMETER side
// 	3.2. bilinear block matching for every block (exact optimization)
// 	(3.2bis. PARAMETER nwarps)
// 4. turbulence compensation (several methods)
//

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#include "fail.c"
#include "xmalloc.c"

static void bilinear_zoom(float *y, int wy, int hy, float *x, int wx, int hx)
{
	// copypaste some standard code here
}

static void laplacian_emphasis_inplace(float *x, int w, int h, float sigma)
{
	// backwards heat equation (some quick and easy implementation)
}

static void cgi_flow(float *f, float *a, float *b, int w, int h,
		int side, int nwarps, float epsilon)
{
	// 1. build, or traverse, patches of side "side"
	//
	// 2. for each patch, solve the corresponding optimization problem
	// (thus producing 6 parameters)
	//
	// 3. for each patch, use the bilinear parameters to produce a dense
	// flow field within that patch
}

static void estimate_the_turbulence(float *y, float **x, int w, int h, int n)
{
}

static void decimate_back(float *y, float *y4, int w, int h)
{
}

void frakes_monaco_smith(float *y, float **x, int w, int h, int n,
		float sigma, int side, int nwarps, float epsilon)
{
	// increase the resolution of each image, and spatially filter them
	float *x4[n];
	for (int i = 0; i < n; i++) {
		x4[i] = xmalloc(4*w * 4*h * sizeof(float));
		bilinear_zoom(x4[i], 4*w, 4*h, x[i], w, h);
		laplacian_emphasis_inplace(x4[i], 4*w, 4*h, sigma);
	}

	// compute the CGI optical flow between consecutive pairs of images
	float *f[n];
	for (int i = 0; i < n-1; i++) {
		f[i] = xmalloc(4*w * 4*h * 2 * sizeof(float));
		cgi_flow(f[i], x4[i], x4[i+1], 4*w, 4*h, side, nwarps, epsilon);
	}

	// estimate the turbulence
	float *y4 = xmalloc(4*w * 4*h * sizeof(float));
	estimate_the_turbulence(y4, x4, w, h, n);
	decimate_back(y, y4, w, h);
}




#ifndef OMIT_MAIN
#include "iio.h"

// auxiliary iterator
static int next(int a, int b, int i)
{
	if (a < b) return i + 1;
	if (a > b) return i - 1;
}

// read a series of images of the same size
static float **read_images(int *w, int *h, int *n, char *s, int a, int b)
{
	if (a == b) fail("I need at least two images!");
	*n = 1 + abs(a - b);
	float **x = xmalloc(*n * sizeof(float));
	int idx = 0;
	*w = *h = 0;
	for (int i = a; i != b; i = next(a,b,i), idx++) {
		int buflen = 0x300, wt, ht;
		char buf[buflen];
		snprintf(buf, buflen, s, i);
		x[idx] = iio_read_image_float(buf, &wt, &ht);
		if ((*w && wt != *w) || (*h && ht != *h))
			fail("input image size mismatch");
		*w = wt;
		*h = ht;
	}
	assert(idx+1 == *n);
	return x;
}

int main(int c, char *v[])
{
	if (c != 8) {
		fprintf(stderr, "usage:\n\t"
		"%s inpat.png first last sigma side nwarps epsilon\n", *v);
	//       0  1         2     3    4     5    6      7
		return EXIT_FAILURE;
	}
	char *inpat = v[1];
	int first_frame = atoi(v[2]);
	int last_frame = atoi(v[3]);
	float sigma = atof(v[4]);
	int   side = atoi(v[5]);
	int   nwarps = atoi(v[6]);
	float epsilon = atof(v[7]);

	int w, h, n;
	float **x = read_images(&n, &w, &h, inpat, first_frame, last_frame);
	assert(n == abs(first_frame - last_frame));

	float *y = xmalloc(w * h * sizeof*y);
	frakes_monaco_smith(y, x, w, h, n, sigma, side, nwarps, epsilon);
	iio_write_image_float("-", y, w, h);

	for (int i = 0; i < n; i++) free(x[i]); free(x); free(y);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
