#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

static void float_sum(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k];
	}
}

static void float_avg(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

static void float_mul(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 1;
		for (int i = 0; i < n; i++)
			y[k] *= x[i][k];
	}
}

static float fnorm(float *x, int n)
{
	switch(n) {
	case 1: return fabs(x[0]);
	case 2: return hypot(x[0], x[1]);
	default: return hypot(x[0], fnorm(x+1, n-1));
	}
}

static void float_min(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float mino = fnorm(x[0], d);
	for (int i = 1; i < n; i++)
	{
		float ni = fnorm(x[i], d);
		if (ni < mino) {
			midx = i;
			mino = ni;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];
}

static void float_max(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float mino = fnorm(x[0], d);
	for (int i = 1; i < n; i++)
	{
		float ni = fnorm(x[i], d);
		if (ni > mino) {
			midx = i;
			mino = ni;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];
}

static float medscore(float *xx, int idx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	float r = 0;
	for (int i = 0; i < n; i++)
		if (i != idx)
		{
			float v[d];
			for (int j = 0; j < d; j++)
				v[j] = x[idx][j] - x[i][j];
			r += fnorm(v, d);
		}
	return r;
}

static void float_med(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float misc = medscore(xx, 0, d, n);
	for (int i = 1; i < n; i++)
	{
		float si = medscore(xx, i, d, n);
		if (si < misc) {
			midx = i;
			misc = si;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];

}

static float float_mod_1d(float *x, int n)
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

//static float float_mod_1dp(float *x, int n, float mi, float ma, int nb)
//{
//	float h[nb];
//	for (int i = 0; i < nb; i++)
//		h[i] = 0;
//	for (int i = 0; i < n; i++)
//	{
//		int xi = x[i];
//		if (xi < 0) fail("negative xi=%g", x[i]);//xi = 0;
//		if (xi > nb-1) fail("large xi=%g", x[i]);//xi = 0xff;
//		h[xi] += 2;
//		if (xi > 0) h[xi-1] += 1;
//		if (xi < nb-1) h[xi+1] += 1;
//	}
//	int midx = nb/2;
//	for (int i = 0; i < nb; i++)
//		if (h[i] > h[midx])
//			midx = i;
//	return midx;
//}

static void float_modc(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int i = 0; i < d; i++) {
		float t[n];
		for (int j = 0; j < n; j++)
			t[j] = x[j][i];
		y[i] = float_mod_1d(t, n);
	}
}

#if 0
static float float_max(float *x, int d, int n)
{
	float r = -INFINIT;
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

static float float_med(float *x, int d, int n)
{
	if (!n) fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return x[0];
	qsort(x, n, sizeof*x, compare_floats);
	return x[n/2];
}

#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

static float float_medv(float *x, int d, int n)
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

static float float_first(float *x, int d, int n)
{
	if (n)
		return x[0];
	else
		fail("empty list of pixel values!");
}

static float float_pick(float *x, int d, int n)
{
	if (n) {
		int i = randombounds(0, n-1);
		return x[i];
	}
	else
		fail("empty list of pixel values!");
}

static bool isgood(float x)
{
	return isfinite(x);
}
#endif

int main(int c, char *v[])
{
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg|mul|med] [v1 ...] > out\n", *v);
		//          0  1                          2  3
		return EXIT_FAILURE;
	}
	int n = c - 2;
	char *operation_name = v[1];
	void (*f)(float*,float*,int,int) = NULL;
	if (0 == strcmp(operation_name, "sum"))   f = float_sum;
	if (0 == strcmp(operation_name, "mul"))   f = float_mul;
	if (0 == strcmp(operation_name, "prod"))  f = float_mul;
	if (0 == strcmp(operation_name, "avg"))   f = float_avg;
	if (0 == strcmp(operation_name, "min"))   f = float_min;
	if (0 == strcmp(operation_name, "max"))   f = float_max;
	if (0 == strcmp(operation_name, "med"))   f = float_med;
	if (0 == strcmp(operation_name, "medi"))   f = float_med;
	if (0 == strcmp(operation_name, "modc"))   f = float_modc;
	//if (0 == strcmp(operation_name, "medv"))   f = float_medv;
	//if (0 == strcmp(operation_name, "rnd"))   f = float_pick;
	//if (0 == strcmp(operation_name, "first")) f = float_first;
	if (!f) fail("unrecognized operation \"%s\"", operation_name);
	float *x[n];
	int w[n], h[n], pd[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_float_vec(v[i+2], w + i, h + i, pd + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h || pd[i] != *pd)
			fail("%dth image sizes mismatch\n", i);
	}
	float (*y) = xmalloc(*w * *h * *pd * sizeof*y);
	for (int i = 0; i < *w * *h; i++)
	{
		float tmp[n][*pd];
		for (int j = 0; j < n; j++)
			for (int k = 0; k < *pd; k++)
				tmp[j][k] = x[j][i**pd+k];
		f(y + i**pd, tmp[0], *pd, n);
	}
	iio_save_image_float_vec("-", y, *w, *h, *pd);
	return EXIT_SUCCESS;
}
