// like "veco", but using doubles
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

static double double_sum(double *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return r;
}

static double double_avg(double *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return n?r/n:r;
}

static double double_mul(double *x, int n)
{
	double r = 1;
	for (int i = 0; i < n; i++)
		r *= x[i];
	return r;
}

static double double_min(double *x, int n)
{
	double r = INFINITY;
	for (int i = 0; i < n; i++)
		if (x[i] < r)
			r = x[i];
	return r;
}

static double double_max(double *x, int n)
{
	double r = -INFINITY;
	for (int i = 1; i < n; i++)
		if (x[i] > r)
			r = x[i];
	return r;
}

int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (*da > *db) - (*da < *db);
}

static double double_med(double *x, int n)
{
	if (!n) return NAN;//fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return x[0];
	qsort(x, n, sizeof*x, compare_doubles);
	return x[n/2];
}

static double double_cnt(double *x, int n)
{
	return n;
}

static double double_mod(double *x, int n)
{
	double h[0x100];
	for (int i = 0; i < 0x100; i++)
		h[i] = 0;
	for (int i = 0; i < n; i++)
	{
		int xi = x[i];
		if (xi < 0) fail("negative xi=%g", x[i]);//xi = 0;
		if (xi > 0xff) fail("large xi=%g", x[i]);//xi = 0xff;
		h[xi] += 2;
		if (xi > 0) h[xi-1] += 1;
		if (xi < 0xff) h[xi+1] += 1;
	}
	int mi = 0x80;
	for (int i = 0; i < 0x100; i++)
		if (h[i] > h[mi])
			mi = i;
	return mi;
}

#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

static double double_medv(double *x, int n)
{
	if (!n) fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return (x[0] + x[1])/2;
	qsort(x, n, sizeof*x, compare_doubles);
	if (EVENP(n))
		return (x[n/2] + x[1+n/2])/2;
	else
		return x[n/2];
}

static double double_percentile(double *x, int n)
{
	static double percentile = 50;
	if (n == -1) { percentile = fmax(0, fmin(*x, 100)); return 0; }
	if (n == 0) return NAN;
	if (n == 1) return x[0];
	qsort(x, n, sizeof*x, compare_doubles);
	int i = round(percentile * (n - 1) / 100.0);
	//fprintf(stderr, "n=%d, percentile=%g, i=%d\n", n, percentile, i);
	assert(i >= 0);
	assert(i < n);
	return x[i];
}

static double double_first(double *x, int n)
{
	if (n)
		return x[0];
	else
		fail("empty list of pixel values!");
}

static double double_pick(double *x, int n)
{
	if (n) {
		int i = randombounds(0, n-1);
		return x[i];
	}
	else
		fail("empty list of pixel values!");
}

typedef bool (*isgood_t)(double);

static bool isgood_finite(double x)
{
	return isfinite(x);
}

static bool isgood_always(double x)
{
	return true;
}

static bool isgood_numeric(double x)
{
	return !isnan(x);
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int gpar = atoi(pick_option(&c, &v, "g", "1"));
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg|mul|med] [v1 ...] > out\n", *v);
		//          0  1                          2  3
		return EXIT_FAILURE;
	}
	int n = c - 2;
	char *operation_name = v[1];
	double (*f)(double *,int) = NULL;
	if (0 == strcmp(operation_name, "sum"))   f = double_sum;
	if (0 == strcmp(operation_name, "mul"))   f = double_mul;
	if (0 == strcmp(operation_name, "prod"))  f = double_mul;
	if (0 == strcmp(operation_name, "avg"))   f = double_avg;
	if (0 == strcmp(operation_name, "min"))   f = double_min;
	if (0 == strcmp(operation_name, "max"))   f = double_max;
	if (0 == strcmp(operation_name, "med"))   f = double_med;
	if (0 == strcmp(operation_name, "mod"))   f = double_mod;
	if (0 == strcmp(operation_name, "cnt"))   f = double_cnt;
	if (0 == strcmp(operation_name, "medi"))   f = double_med;
	if (0 == strcmp(operation_name, "medv"))   f = double_medv;
	if (0 == strcmp(operation_name, "rnd"))   f = double_pick;
	if (0 == strcmp(operation_name, "first")) f = double_first;
	if (*operation_name == 'q') {
		double p = atof(1 + operation_name);
		f = double_percentile;
		f(&p, -1);
	}
	if (!f) fail("unrecognized operation \"%s\"", operation_name);
	bool (*isgood)(double) = NULL;
	if (0 == gpar) isgood = isgood_finite;
	if (1 == gpar) isgood = isgood_numeric;
	if (2 == gpar) isgood = isgood_always;
	if (!isgood) fail("unrecognized goodness %d", gpar);
	double *x[n];
	int w[n], h[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_double(v[i+2], w + i, h + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h)
			fail("%dth image size mismatch\n", i);
	}
	double (*y) = xmalloc(*w * *h * sizeof*y);
	for (int i = 0; i < *w * *h; i++)
	{
		double tmp[n];
		int ngood = 0;
		for (int j = 0; j < n; j++)
			if (isgood(x[j][i]))
				tmp[ngood++] = x[j][i];
		y[i] = f(tmp, ngood);
	}
	iio_write_image_double("-", y, *w, *h);
	return EXIT_SUCCESS;
}
