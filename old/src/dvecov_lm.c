#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "pickopt.c"

static void acc_euc(double *a, double *x, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = hypot(a[i], x[i]);
}

static void acc_sum(double *a, double *x, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = a[i] + x[i];
}

static void acc_mul(double *a, double *x, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = a[i] * x[i];
}

static void acc_avg(double *a, double *x, int n)
{
	acc_sum(a, x, n); // post-processed later
}

static void acc_min(double *a, double *x, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = fmin(a[i], x[i]);
}

static void acc_max(double *a, double *x, int n)
{
	for (int i = 0; i < n; i++)
		a[i] = fmax(a[i], x[i]);
}

static double euclidean_norm(double *x, int n)
{
	return n ? hypot(*x, euclidean_norm(x+1, n-1)) : 0;
}

static void acc_Min(double *a, double *x, int n)
{
	double na = euclidean_norm(a, n);
	double nx = euclidean_norm(x, n);
	for (int i = 0; i < n; i++)
		a[i] = na < nx ? a[i] : x[i];
}

static void acc_Max(double *a, double *x, int n)
{
	double na = euclidean_norm(a, n);
	double nx = euclidean_norm(x, n);
	for (int i = 0; i < n; i++)
		a[i] = na > nx ? a[i] : x[i];
}

static bool isgood(double *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isfinite(x[i]))
			return false;
	return true;
}

int main(int c, char *v[])
{
	// process input arguments
	char *filename_out = pick_option(&c, &v, "o", "-");
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg} [v1 ...] [-o out]\n", *v);
		//          0  1                 2
		return 1;
	}
	int n = c - 2;
	char *operation_name = v[1];
	char **filename_in = v+2;

	// identify operator
	void (*f)(double*,double*,int) = NULL;
	if (0 == strcmp(operation_name, "euc"))   f = acc_euc;
	if (0 == strcmp(operation_name, "sum"))   f = acc_sum;
	if (0 == strcmp(operation_name, "avg"))   f = acc_avg;
	if (0 == strcmp(operation_name, "min"))   f = acc_min;
	if (0 == strcmp(operation_name, "max"))   f = acc_max;
	if (0 == strcmp(operation_name, "Min"))   f = acc_Min;
	if (0 == strcmp(operation_name, "Max"))   f = acc_Max;
	if (0 == strcmp(operation_name, "mul"))   f = acc_mul;
	if (!f) fail("unrecognized operation \"%s\"", operation_name);

	// read first image
	int w, h, pd;
	double *x = iio_read_image_double_vec(filename_in[0], &w, &h, &pd);

	// accumulate the other images
	for (int i = 1; i < n; i++)
	{
		int W, H, PD;
		double *y =iio_read_image_double_vec(filename_in[i], &W, &H, &PD);
		if (W!=w)fail("%dth image width mismatch %d != %d\n",i,w,W);
		if (H!=h)fail("%dth image width mismatch %d != %d\n",i,h,H);
		if (PD!=pd)fail("%dth image depth mismatch %d != %d\n",i,pd,PD);
		for (int j = 0; j < w*h; j++)
			f(x + pd * j, y + pd * j, pd);
		free(y);
	}
	if (f == acc_avg)
		for (int i = 0; i < w*h*pd; i++)
			x[i] /= n;

	// save accumulated image
	iio_write_image_double_vec(filename_out, x, w, h, pd);


	// cleanup and exit
	free(x);
	return 0;
}
