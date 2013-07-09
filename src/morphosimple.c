#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new) {
		fprintf(stderr, "ERROR: out of memory when requesting "
			       "%zu bytes\n", size);
		abort();
	}
	return new;
}

typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

// extrapolate by nan
static float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return NAN;
	return x[i + j*w];
}

// extrapolate by nearest value
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i + j*w];
}


// a structuring element is a list E of integers
// E[0] = number of pixels
// E[1] = 0 (flags, not used yet)
// (E[2], E[3]) = position of the center
// (E[4], E[5]) = first pixel
// (E[6], E[7]) = second pixel
// ...


void morsi_erosion(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmin(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

void morsi_dilation(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = -INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmax(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

float median(float *a, int n)
{
	if (n < 1) return NAN;
	if (n == 1) return *a;
	if (n == 2) return (a[0] + a[1])/2;
	qsort(a, n, sizeof*a, compare_floats);
	if (0 == n%2)
		return (a[n/2]+a[1+n/2])/2;
	else
		return a[n/2];
}

void morsi_median(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a[e[0]];
		int cx = 0;
		for (int k = 0; k < e[0]; k++)
		{
			float v = p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]);
			if (isfinite(v))
				a[cx++] = v;
		}
		y[j*w+i] = median(a, cx);
	}
}

void morsi_opening(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	morsi_dilation(y, t, w, h, e);
	free(t);
}

void morsi_closing(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	morsi_dilation(y, t, w, h, e);
	free(t);
}

void morsi_gradient(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_erosion(a, x, w, h, e);
	morsi_dilation(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = b[i] - a[i];
	free(a);
	free(b);
}

#define MORSI_TEST_MAIN

#ifdef MORSI_TEST_MAIN
#include <string.h>
#include "iio.h"
int main(int c, char **v)
{
	// data for available structuring elements
	int cross[] = {5,0,  0,0, -1,0, 0,0, 1,0, 0,-1, 0,1 };
	int square[] = {9,0, 0,0, -1,-1,-1,0,-1,1, 0,-1,0,0,0,1, 1,-1,1,0,1,1};

	// process input arguments
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
				"%s element operation [in [out]]\n", *v);
		//                0 1       2          3   4
		return 1;
	}
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";
	int *structuring_element = NULL;
	if (0 == strcmp(v[1], "cross" )) structuring_element = cross;
	if (0 == strcmp(v[1], "square")) structuring_element = square;
	if (!structuring_element) {
		fprintf(stderr, "elements = cross, square ...\n");
		return 1;
	}
	void (*operation)(float*,float*,int,int,int*) = NULL;
	if (0 == strcmp(v[2], "erosion"  )) operation = morsi_erosion;
	if (0 == strcmp(v[2], "dilation" )) operation = morsi_dilation;
	if (0 == strcmp(v[2], "median"   )) operation = morsi_median;
	if (0 == strcmp(v[2], "opening"  )) operation = morsi_opening;
	if (0 == strcmp(v[2], "closing"  )) operation = morsi_closing;
	if (0 == strcmp(v[2], "gradient" )) operation = morsi_gradient;
	if (!operation) {
		fprintf(stderr, "operations = erosion, dilation, opening...\n");
		return 1;
	}

	// prepare input and output images
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *y = malloc(w*h*sizeof*y);

	// compute
	operation(y, x, w, h, structuring_element);

	// save result
	iio_save_image_float(filename_out, y, w, h);

	// cleanup
	free(x);
	free(y);
	return 0;
}
#endif//MORSI_TEST_MAIN
