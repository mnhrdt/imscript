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

static float float_sum(float *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return r;
}

static float float_avg(float *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return n?r/n:r;
}

static float float_mul(float *x, int n)
{
	double r = 1;
	for (int i = 0; i < n; i++)
		r *= x[i];
	return r;
}

static float float_min(float *x, int n)
{
	float r = INFINITY;
	for (int i = 0; i < n; i++)
		if (x[i] < r)
			r = x[i];
	return r;
}

static float float_max(float *x, int n)
{
	float r = -INFINITY;
	for (int i = 1; i < n; i++)
		if (x[i] > r)
			r = x[i];
	return r;
}

int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static float float_med(float *x, int n)
{
	if (!n) return NAN;//fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return x[0];
	qsort(x, n, sizeof*x, compare_floats);
	return x[n/2];
}

static float float_cnt(float *x, int n)
{
	return n;
}

static float float_mod(float *x, int n)
{
	float h[0x100];
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

static float float_medv(float *x, int n)
{
	if (!n) fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return (x[0] + x[1])/2;
	qsort(x, n, sizeof*x, compare_floats);
	if (EVENP(n))
		return (x[n/2] + x[1+n/2])/2;
	else
		return x[n/2];
}

static float float_percentile(float *x, int n)
{
	static float percentile = 50;
	if (n == -1) { percentile = fmax(0, fmin(*x, 100)); return 0; }
	if (n == 0) return NAN;
	if (n == 1) return x[0];
	qsort(x, n, sizeof*x, compare_floats);
	int i = round(percentile * (n - 1) / 100.0);
	//fprintf(stderr, "n=%d, percentile=%g, i=%d\n", n, percentile, i);
	assert(i >= 0);
	assert(i < n);
	return x[i];
}

static float float_first(float *x, int n)
{
	if (n)
		return x[0];
	else
		fail("empty list of pixel values!");
}

static float float_pick(float *x, int n)
{
	if (n) {
		int i = randombounds(0, n-1);
		return x[i];
	}
	else
		fail("empty list of pixel values!");
}

typedef bool (*isgood_t)(float);

static bool isgood_finite(float x)
{
	return isfinite(x);
}

static bool isgood_always(float x)
{
	return true;
}

static bool isgood_numeric(float x)
{
	return !isnan(x);
}

#include "pickopt.c"

int main_veco(int c, char *v[])
{
	int gpar = atoi(pick_option(&c, &v, "g", "1"));
	char *filename_out = pick_option(&c, &v, "o", "-");
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg|mul|med] [v1 ...] > out\n", *v);
		//          0  1                          2  3
		return EXIT_FAILURE;
	}
	int n = c - 2;
	char *operation_name = v[1];
	float (*f)(float *,int) = NULL;
	if (0 == strcmp(operation_name, "sum"))   f = float_sum;
	if (0 == strcmp(operation_name, "mul"))   f = float_mul;
	if (0 == strcmp(operation_name, "prod"))  f = float_mul;
	if (0 == strcmp(operation_name, "avg"))   f = float_avg;
	if (0 == strcmp(operation_name, "min"))   f = float_min;
	if (0 == strcmp(operation_name, "max"))   f = float_max;
	if (0 == strcmp(operation_name, "med"))   f = float_med;
	if (0 == strcmp(operation_name, "mod"))   f = float_mod;
	if (0 == strcmp(operation_name, "cnt"))   f = float_cnt;
	if (0 == strcmp(operation_name, "medi"))   f = float_med;
	if (0 == strcmp(operation_name, "medv"))   f = float_medv;
	if (0 == strcmp(operation_name, "rnd"))   f = float_pick;
	if (0 == strcmp(operation_name, "first")) f = float_first;
	if (*operation_name == 'q') {
		float p = atof(1 + operation_name);
		f = float_percentile;
		f(&p, -1);
	}
	if (!f) fail("unrecognized operation \"%s\"", operation_name);
	bool (*isgood)(float) = NULL;
	if (0 == gpar) isgood = isgood_finite;
	if (1 == gpar) isgood = isgood_numeric;
	if (2 == gpar) isgood = isgood_always;
	if (!isgood) fail("unrecognized goodness %d", gpar);
	float *x[n];
	int w[n], h[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_float(v[i+2], w + i, h + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h)
			fail("%dth image size mismatch\n", i);
	}
	float (*y) = xmalloc(*w * *h * sizeof*y);
	for (int i = 0; i < *w * *h; i++)
	{
		float tmp[n];
		int ngood = 0;
		for (int j = 0; j < n; j++)
			if (isgood(x[j][i]))
				tmp[ngood++] = x[j][i];
		y[i] = f(tmp, ngood);
	}
	iio_write_image_float(filename_out, y, *w, *h);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_veco(c, v); }
#endif
