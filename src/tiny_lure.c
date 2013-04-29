#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>


struct floatint { float f; int i; } ;
static int compare_struct_floatint(const void *aa, const void *bb)
{
	const struct floatint *a = (const struct floatint *)aa;
	const struct floatint *b = (const struct floatint *)bb;
	return (a->f - b->f) - (a->f < b->f);
}

static int winpatch_length(int pd, float W)
{
	int wi = W;
	int wis = 2*wi+1;
	int r = pd*wis*wis;
	return r;
}

#include "getpixel.c"

static void getpatch(float *p, float *x, int w, int h, int pd,
                                         int i, int j, float W)
{
	int rad = W;
	int cx = 0;
	for (int di = -rad; di <= rad; di++)
	for (int dj = -rad; dj <= rad; dj++)
	for (int l = 0; l < pd; l++)
		p[cx++] = getsample_0(x, w, h, pd, i+di, j+dj, l);
	assert(cx == winpatch_length(pd, W));
}

static float patchdistance(float *p, float *q, int n)
{
	float r = 0;
	for (int i = 0; i < n; i++)
		r = hypot(r, p[i] - q[i]);
	return r;
}

static float wpatch_distance(float *x, float *y, int w, int h, int pd,
                                                 int i, int j, float W)
{
	int nW = winpatch_length(pd, W);
	float px[nW], py[nW];
	getpatch(px, x, w, h, pd, i, j, W);
	getpatch(py, y, w, h, pd, i, j, W);
	return patchdistance(px, py, nW);
}

#define OMIT_BLUR_MAIN
#include "blur.c"

void silly_lucky_region(float **y, float **x, int n, int w, int h, int pd,
		float W, float S)
{
	float *mx = malloc(w*h*pd*sizeof*mx);
	for (int i = 0; i < w*h*pd; i++)
	{
		float m = 0;
		for (int j = 0; j < n; j++)
			m += x[j][i];
		mx[i] = m/n;
	}

	float *bx[n], blurparams[1] = {S};
	for (int i = 0; i < n; i++)
	{
		bx[i] = malloc(w * h * pd * sizeof*bx[0]);
		blur_2d(bx[i], x[i], w, h, pd, "gaussian", blurparams, 1);
	}

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		struct floatint v[n];
		for (int k = 0; k < n; k++)
		{
			v[k].i = k;
			v[k].f = wpatch_distance(mx, bx[k], w, h, pd, i, j, W);
		}
		qsort(v, n, sizeof*v, compare_struct_floatint);
		int idx = j*w + i;
		for (int k = 0; k < n; k++)
			for (int l = 0; l < pd; l++)
				y[k][idx*pd+l] = x[v[k].i][idx*pd+l];
	}

	free(mx);
	for (int i = 0; i < n; i++)
		free(bx[i]);
}

#include "iio.h"

#include "fail.c"

int main(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s winsize ker inpat outpat first last\n", *v);
		//        0 1       2   3     4      5     6
		return 1;
	}
	float winsize = atof(v[1]);
	float kersigma = atof(v[2]);
	char *filepattern_in = v[3];
	char *filepattern_out = v[4];
	int idx_first = atoi(v[5]);
	int idx_last = atoi(v[6]);

	int n = idx_last - idx_first + 1;
	if (n < 3 || n > 1000) fail("bad n = %d\n", n);

	int w, h, pd = 0;
	float *x[n], *y[n];
	for (int i = 0; i < n; i++)
	{
		int idx = idx_first + i;
		char filename_in[FILENAME_MAX];
		snprintf(filename_in, FILENAME_MAX, filepattern_in, idx);
		int ww, hh, ppdd;
		x[i] = iio_read_image_float_vec(filename_in, &ww, &hh, &ppdd);
		if (pd != 0) {
			if (w != ww || h != hh || pd != ppdd)
				fail("input images size mismatch");
		} else {
			w = ww;
			h = hh;
			pd = ppdd;
		}
	}

	for (int i = 0; i < n; i++)
		y[i] = malloc(w * h * pd * sizeof*y[0]);

	silly_lucky_region(y, x, n, w, h, pd, winsize, kersigma);

	for (int i = 0; i < n; i++)
	{
		int idx = idx_first + i;
		char filename_out[FILENAME_MAX];
		snprintf(filename_out, FILENAME_MAX, filepattern_out, idx);
		iio_save_image_float_vec(filename_out, y[i], w, h, pd);
	}

	return 0;
}
