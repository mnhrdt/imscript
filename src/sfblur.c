#include <assert.h>
#include "xmalloc.c"

static float getcomponent_0(float *x, int n, int i)
{
	return (i < 0 || i >= n) ? 0 : x[i];
}

static float getcomponent_1(float *x, int n, int i)
{
	if (i < 0) i = 0;
	if (i >= n) i = n - 1;
	return x[i];
}

// like n%p, but works for all numbers
static int good_modulus(int n, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(n, -p);

	int r = n % p;
	r = r < 0 ? r + p : r;

	assert(r >= 0);
	assert(r < p);
	return r;
}

// symmetrized and periodized index
static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p)
		r -= 1;
	if (r > p)
		r = 2*p - r;
	assert(r >= 0);
	assert(r < p);
	return r;
}

static float getcomponent_2(float *x, int n, int i)
{
	i = positive_reflex(i, n);
	return x[i];
}

void slow_blur_1d(float *x, int n, float side)
{
	int k = side;  // integer-valued side
	int o = k / 2; // left-offset
	float t[n];    // temporary storage
	for (int i = 0; i < n; i++)
	{
		t[i] = 0;
		for (int j = 0; j < k; j++)
			t[i] += getcomponent_1(x, n, i + j - o);
	}
	for (int i = 0; i < n; i++)
		x[i] = t[i]/k;
}

void fast_blur_1d(float *x, int n, float side)
{
	int k = side;  // integer-valued side
	int o = k / 2; // left-offset
	float t[n];    // temporary storage
	// compute the first element
	t[0] = o*x[0];
	for (int i = 1; i < o; i++)
		t[0] += x[i];
	for (int i = 1; i < n; i++)
	{
		float treu = getcomponent_2(x, n, i - o);
		float posa = getcomponent_2(x, n, i - o + k);
		t[i] = t[i-1] - treu + posa;
	}
	for (int i = 0; i < n; i++)
		x[i] = t[i]/k;
}

void superfast_blur_gray(float *x, int w, int h, float side)
{
	// blur in the horizontal direction
	for (int j = 0; j < h; j++)
		slow_blur_1d(x + j*w, w, side);
		//fast_blur_1d(x + j*w, w, side);

	// blur in the vertical direction
	// (TODO)
}

void superfast_blur_split(float *x, int w, int h, int pd, float side)
{
	for (int i = 0; i < pd; i++)
		superfast_blur_gray(x + w*h*i, w, h, side);
}

#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s side [in [out]]\n", *v);
		//                0 1      2   3
		return 1;
	}
	float side = atof(v[1]);
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// perform the computation
	superfast_blur_split(x, w, h, pd, side);

	// save output image
	iio_write_image_float_split(filename_out, x, w, h, pd);

	// cleanup and exit
	free(x);
	return 0;
}
