#include <stdbool.h>
#include <stdio.h>
#include <math.h>



#include "fail.c"
#include "xmalloc.c"
#include "smapa.h"


static int get_nvals(float *v, float *wv, float *x, int w, int h, int i, int j)
{
	int r = 0, n[8][2] = {
		{1,0}, {0,1}, {-1,0}, {0,-1},
		{1,1}, {-1,1}, {-1,-1}, {1,-1}
	};
	for (int p = 0; p < 4; p++)
	{
		int ii = i + n[p][0];
		int jj = j + n[p][1];
		if (ii >= 0 && jj >= 0 && ii < w && jj < h)
		{
			v[r] = x[w*jj+ii];
			if (wv)
				wv[r] = hypot(ii,jj);
			r += 1;
		}
	}
	return r;
}

static void get_minmax(float *min, float *max, float *x, int n)
{
	*min = *max = x[0];
	for (int i = 1; i < n; i++)
	{
		if (x[i] < *min) *min = x[i];
		if (x[i] > *max) *max = x[i];
	}
}


SMART_PARAMETER(AMLE_NITER,100)

// the input mask is coded by nans
void amle(float *y, float *x, int w, int h)
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
			y[idx] = 0;
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
			int idx = j*w + i;
			float v[0x100], wv[0x100], mi, ma;
			int nv = get_nvals(v, NULL, y, w, h, i, j);
			get_minmax(&mi, &ma, v, nv);
			float newy = (mi + ma)/2;
			if (fabs(newy - y[idx]) > actumax)
				actumax = fabs(newy - y[idx]);
			y[idx] = newy;
		}
		if (0 == iter % 10)
			fprintf(stderr, "iter %d, actu = %g\n", iter, actumax);
	}

	free(mask);
}

#include "iio.h"

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

	iio_save_image_float(filename_out, out, *w, *h);

	return 0;
}
