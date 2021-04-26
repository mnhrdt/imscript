#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "iio.h"

static float weighted_average(float *x, float *w, int n)
{
	double m = 0; // weighted sum
	double r = 0; // sum of weights
	for (int i = 0; i < n; i++)
	{
		m += w[i] * x[i];
		r += w[i];
	}
	return m / r;
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel_1(float *I, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return I[i+j*w];
}


// mediator with weights
void mediatorw(
		float *y,   // output image
		float *x,   // input image
		int w,      // image width
		int h,      // image height
		int r,      // filter radius
		int n       // number of reweighting iterartions
		)
{
	// kernel side
	int s = (2*r + 1);
	fprintf(stderr, "s = %d, s^2 = %d\n", s, s*s);

	// initialize all weigths to 1
	float *W = malloc(w*h * s*s * sizeof*W);
	for (int i = 0; i < w*h * s*s; i++)
		W[i] = 1;

	// iterate n times
	for (int iter = 0; iter < n; iter++)
	{
		// compute weighted average
#ifdef _OPENMP
#pragma omp parallel for
#endif//_OPENMP
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			int cx = 0;   // counter
			float t[s*s]; // patch of x around (i,j)
			for (int l = -r; l <= r; l++)
			for (int k = -r; k <= r; k++)
				t[cx++] = getpixel_1(x, w, h, i+k, j+l);
			assert(cx == s*s);
			y[j*w+i] = weighted_average(t, W+s*s*(j*w+i), s*s);
		}

		// update weights
#ifdef _OPENMP
#pragma omp parallel for
#endif//_OPENMP
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			int cx = 0;   // counter
			float t[s*s]; // patch of x around (i,j)
			for (int l = -r; l <= r; l++)
			for (int k = -r; k <= r; k++)
				t[cx++] = getpixel_1(x, w, h, i+k, j+l);
			assert(cx == s*s);
			int ij = j*w+i;
			for (int k = 0; k < s*s; k++)
				W[s*s*ij+k] = 1 / (1e-6 + fabs(t[k] - y[ij]));
		}

	}

	// cleanup and exit
	free(W);
}

void mediatorw_separable(float *y, float *x, int w, int h, int pd, int r, int n)
{
	for (int k = 0; k < pd; k++)
		mediatorw(y + k*w*h, x + k*w*h, w, h, r, n);
}

#include <stdio.h>
int main(int c, char *v[])
{
	if ( c < 3 ) {
		fprintf(stderr, "usage:\n\t%s r n [in [out]]\n", *v);
		//                          0 1 2  3   4
		return 0;
	}
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";
	int r = atoi(v[1]); // radius
	int n = atoi(v[2]); // number of iterations

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = malloc(w * h * pd * sizeof*y);

	mediatorw_separable(y, x, w, h, pd, r, n);

	iio_write_image_float_split(filename_out, y, w, h, pd);

	return 0;
}
