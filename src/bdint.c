// lowest neighbor interpolation



#include <assert.h>
#include <stdlib.h>
#include <math.h>


static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

#define xmalloc malloc
#include "abstract_dsf.c"

void bdint(float *x, int w, int h)
{
	// create the dsf
	int *rep = xmalloc(w*h*sizeof*rep);
	adsf_begin(rep,w*h);

	// remove from dsf pixels with known values in input
	for (int i = 0; i < w*h; i++)
		if (!isnan(x[i]))
			rep[i] = -1;

	// join neighboring NANs
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (isnan(x[p0]) && isnan(x[p1]))
			adsf_union(rep, w*h, p0, p1);
		if (isnan(x[p0]) && isnan(x[p2]))
			adsf_union(rep, w*h, p0, p2);
	}

	// canonicalize dsf
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			rep[i] = adsf_find(rep, w*h, i);

	// prepare table of optima
	float *opt = xmalloc(w*h*sizeof*opt);
	for (int i = 0; i < w*h; i++)
		opt[i] = INFINITY;

	// for each valued point that has a neighboring nan, update the optimum
	int n[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int k = 0; k < 4; k++)
	{
		int ii = i + n[k][0];
		int jj = j + n[k][1];
		int idx0 = j  * w + i;
		int idx1 = jj * w + ii;
		if (insideP(w, h, i, j) && insideP(w, h, ii, jj) &&
				!isnan(x[idx0]) && isnan(x[idx1]))
			opt[rep[idx1]] = fmin(opt[rep[idx1]], x[idx0]);
	}

	// fill-in the computed optima
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ij = j  * w + i;
		if (isnan(x[ij]))
			x[ij] = opt[rep[ij]];
	}

	//cleanup
	free(opt);
	free(rep);
}


#ifndef OMIT_BDINT_MAIN
#define USE_BDINT_MAIN
#endif

#ifdef USE_BDINT_MAIN
#include <stdio.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s in.tiff out.tiff]\n", *v);
		//                          0 1       2
		return 1;
	}
	char *filename_in   = v[1];
	char *filename_out  = v[2];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	bdint(x, w, h);

	iio_save_image_float(filename_out, x, w, h);

	return 0;
}
#endif//USE_BDINT_MAIN
