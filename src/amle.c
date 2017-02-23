#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>


#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "smapa.h"

SMART_PARAMETER(AMLE_NN,4)

static int get_nvals(float *v, float *wv2, float *x, int w, int h, int i, int j)
{
	int r = 0, n[][3] = { // {x, y, x*x+y*y}
		{+1,0,1}, {0,+1,1}, {-1,0,1}, {0,-1,1}, // 4-connexity
		{+1,+1,2}, {-1,-1,2}, // 6-connexity
		{-1,+1,2}, {+1,-1,2}, // 8-connexity
		{+2,+1,5}, {+1,+2,5}, {+2,-1,5}, {+1,-2,5},
		{-2,-1,5}, {-1,-2,5}, {-2,+1,5}, {-1,+2,5}, // 16-neighbors
		{+3,+1,10}, {+1,+3,10}, {+3,-1,10}, {+1,-3,10},
		{-3,-1,10}, {-1,-3,10}, {-3,+1,10}, {-1,+3,10}, // 24-neighbors
		{+3,+2,13}, {+2,+3,13}, {+3,-2,13}, {+2,-3,13},
		{-3,-2,13}, {-2,-3,13}, {-3,+2,13}, {-2,+3,13}, // 32-neighbors
		//{+2,0,4}, {0,+2,4}, {-2,0,4}, {0,-2,4}, (non-primitive)
		//{+4,0,16}, {0,+4,16}, {-4,0,16}, {0,-4,16} (non-primitive)
	};
	int nn = AMLE_NN();
	for (int p = 0; p < nn; p++)
	{
		int ii = i + n[p][0];
		int jj = j + n[p][1];
		if (ii >= 0 && jj >= 0 && ii < w && jj < h)
		{
			v[r] = x[w*jj+ii];
			if (wv2)
				wv2[r] = n[p][2];
			r += 1;
		}
	}
	return r;
}

static void get_minmax(float *min, float *max, float *x, int n)
{
	*min = INFINITY;
	*max = -INFINITY;
	for (int i = 0; i < n; i++)
	{
		if (x[i] < *min) *min = x[i];
		if (x[i] > *max) *max = x[i];
	}
}

static void get_minmax_idx(int *min, int *max, float *x, int n)
{
	*min = *max = 0;
	for (int i = 1; i < n; i++)
	{
		if (x[i] < x[*min]) *min = i;
		if (x[i] > x[*max]) *max = i;
	}
}

static float amle_iteration(float *x, int w, int h, int (*mask)[2], int nmask)
{
	float actus = 0;
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i, min, max;
		float value[0x100], weight[0x100];
		int nv = get_nvals(value, weight, x, w, h, i, j);
		get_minmax_idx(&min, &max, value, nv);
		float a = weight[max];
		float b = weight[min];
		float newx = (a*value[min] + b*value[max]) / (a + b);
		actus += fabs(x[idx] - newx);
		//if (fabs(x[idx]-newx) > actumax)
		//	actumax = fabs(x[idx]-newx);
		x[idx] = newx;
	}
	return actus;
}

static void amle_init(float *tinf, float *tsup, float *x, int w, int h)
{
	float min, max;
	get_minmax(&min, &max, x, w*h);
	for (int i = 0; i < w*h; i++)
	{
		if (isnan(x[i])) {
			tinf[i] = min;
			tsup[i] = max;
		} else {
			tinf[i] = x[i];
			tsup[i] = x[i];
		}
	}
	char *filename_init = getenv("AMLE_INIT");
	if (filename_init) {
		int ww, hh;
		float *ini = iio_read_image_float(filename_init, &ww, &hh);
		if (w*h != ww*hh) fail("init size mismatch");
		for (int i = 0; i < w*h; i++)
			if (isnan(x[i]))
				tsup[i] = ini[i];
		free(ini);
	}
}

static float absolute_difference(float *a, float *b, int w, int h,
		int (*mask)[2], int nmask)
{
	float r = 0;
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i;
		float t = fabs(a[idx] - b[idx]);
		if (t > r)
			r = t;
	}
	return r;
}

static float mean_difference(float *a, float *b, int w, int h,
		int (*mask)[2], int nmask)
{
	float r = 0;
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i;
		r += fabs(a[idx] - b[idx]);
	}
	return r/nmask;
}

SMART_PARAMETER(AMLE_TAU,0.25)

static void amle_refine(float *a, float *b, int w, int h,
		int (*mask)[2], int nmask)
{
	float t = AMLE_TAU();
	assert(t >= 0);
	assert(t < 0.5);
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i;
		float newa = (1-t)*a[idx] + t*b[idx];
		float newb = (1-t)*b[idx] + t*a[idx];
		a[idx] = newa;
		b[idx] = newb;
	}
}

//static int randombounds(int a, int b)
//{
//	if (b < a)
//		return randombounds(b, a);
//	if (b == a)
//		return b;
//	return a + rand()%(b - a + 1);
//}

//static void swap(void *a, void *b, size_t s)
//{
//#if 0
//#error "memcpy is way slower!"
//	char t[s];
//	memcpy(t, a, s);
//	memcpy(a, b, s);
//	memcpy(b, t, s);
//#else
//	char *x = a;
//	char *y = b;
//	for (unsigned int i = 0; i < s; i++, x++, y++)
//	{
//		char t = *x;
//		*x = *y;
//		*y = t;
//	}
//#endif
//}

//static void shuffle(void *t, int n, size_t s)
//{
//	char *c = t;
//
//	for (int i = 0; i < n-1; i++)
//		swap(c + s*i, c + s*randombounds(i, n-1), s);
//}

SMART_PARAMETER(AMLE_NITER,100)

void amle(float *y, float *x, int w, int h)
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
		}

	float *tinf = xmalloc(w*h*sizeof(float));
	float *tsup = xmalloc(w*h*sizeof(float));
	amle_init(tinf, tsup, x, w, h);

	for (int iter = 0 ; iter < AMLE_NITER(); iter++)
	{
		float actus_inf = amle_iteration(tinf, w, h, mask, nmask);
		float actus_sup = amle_iteration(tsup, w, h, mask, nmask);
		float e = absolute_difference(tinf, tsup, w, h, mask, nmask);
		float ea = mean_difference(tinf, tsup, w, h, mask, nmask);

		if (0 == iter % 10)
			fprintf(stderr,
				"iter %d, e = {%g %g}, eactus = {%g %g}\n",
				iter, e, ea, actus_inf, actus_sup);

		//if (0 == iter % 33)
		//	shuffle(mask, nmask, sizeof*mask);

		//if (0 == iter % 10)
		//	amle_refine(tinf, tsup, w, h, mask, nmask);
	}

	free(mask);

	for (int i = 0; i < w*h; i++)
		//y[i] = (tinf[i] + tsup[i])/2;
		y[i] = tsup[i];
}

// the input mask is coded by nans
#if 0
void amle_old(float *y, float *x, int w, int h)
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	fprintf(stderr, "nmask = %d\n", nmask);
	int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;
		if (isnan(x[idx])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
			y[idx] = AMLE_INIT();
		} else
			y[idx] = x[idx];
	}

	for (int iter = 0 ; iter < AMLE_NITER(); iter++)
	{
		float actumax = 0;
		for (int p = 0; p < nmask; p++)
		{
			int i = mask[p][0];
			int j = mask[p][1];
			int idx = j*w + i, mi, ma;
			float v[0x100], wv[0x100];
			int nv = get_nvals(v, wv, y, w, h, i, j);
			get_minmax_idx(&mi, &ma, v, nv);
			//float newy = (mi + ma)/2;
			float newy = wv[ma]*v[mi] + wv[mi]*v[ma];
			newy /= wv[mi] + wv[ma];
			if (fabs(newy - y[idx]) > actumax)
				actumax = fabs(newy - y[idx]);
			y[idx] = newy;
		}
		if (0 == iter % 10)
			fprintf(stderr, "iter %d, actu = %g\n", iter, actumax);
	}

	free(mask);
}
#endif


int main(int c, char *v[])
{
	if (c > 4 || (v[1] && v[1][0]=='-' && v[1][1]=='h')) {
		fprintf(stderr, "usage:\n\t%s [in [mask [out]]]\n", *v);
		//                          0  1   2     3
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_mask = c > 2 ? v[2] : "";
	char *filename_out = c > 3 ? v[3] : "-";

	int w[2], h[2];
	float *in = iio_read_image_float(filename_in, w, h);
	float *mask = NULL;
	if (filename_mask[0]) {
		mask = iio_read_image_float(filename_mask, w+1, h+1);
		if (w[0] != w[1] || h[0] != h[1])
			fail("image and mask file size mismatch");
	}
	float *out = xmalloc(*w**h*sizeof*out);

	if (mask)
		for (int i = 0; i < *w**h; i++)
			if (mask[i] > 0)
				in[i] = NAN;

	amle(out, in, *w, *h);

	iio_write_image_float(filename_out, out, *w, *h);

	return 0;
}
