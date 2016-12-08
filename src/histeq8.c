#include <math.h>
#include <stdio.h>

static double transport_distance(double *a, double *b, int n)
{
	// accumulate
	double aa[n], bb[n];
	aa[0] = a[0];
	bb[0] = b[0];
	for (int i = 1; i < n; i++)
	{
		aa[i] = aa[i-1] + a[i];
		bb[i] = bb[i-1] + b[i];
	}

	// normalize
	for (int i = 0; i < n; i++)
	{
		aa[i] /= aa[n-1];
		bb[i] /= bb[n-1];
	}

	// compute L1 norm of the difference
	double r = 0;
	for (int i = 0; i < n; i++)
		r += fabs(aa[i] - bb[i]);

	return r;
}

static void accumulate_histogram(long double h[256])
{
	for (int i = 1; i < 256; i++)
		h[i] += h[i-1];
}

static void fill_histogram(long double h[256], float *x, int n)
{
	for (int i = 0; i < 256; i++)
		h[i] = 0;
	for (int i = 0; i < n; i++)
	{
		int idx = round(x[i]);
		if (idx < 0) idx = 0;
		if (idx > 255) idx = 255;
		h[idx] += 1;
	}
}


static void equalize_inplace(float *x, long double h[256], int n)
{
	for (int i = 0; i < n; i++)
	{
		int idx = round(x[i]);
		if (idx < 0) idx = 0;
		if (idx > 255) idx = 255;
		x[i] = h[idx]/n;
	}
}

#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                         0   1   2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	long double H[256];
	fill_histogram(H, x, w*h);
	accumulate_histogram(H);
	equalize_inplace(x, H, w*h);

	iio_save_image_float(filename_out, x, w, h);

	return 0;
}
