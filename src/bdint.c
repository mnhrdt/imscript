// lowest neighbor interpolation
// optionally: highest-neighbor and average-neighbor
// todo: IDW/shepard of all the neighbors



#include <assert.h>
#include <stdlib.h>
#include <math.h>


static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

#define xmalloc malloc
#include "abstract_dsf.c"

typedef float (accumulator_t)(float, float);

static float accumulate_min(float x, float y)
{
	return fmin(x, y);
}

static float accumulate_max(float x, float y)
{
	return fmax(x, y);
}

static float counted_accumulator(float a, float x, int n)
{
	return (a*n+x)/(n+1);
}

void bdint(float *x, int w, int h, accumulator_t *a)
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

	// canonicalize dsf (after this, the DSF is not changed anymore)
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			rep[i] = adsf_find(rep, w*h, i);

	// prepare table of optima
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i)
			x[i] = -a(-INFINITY,INFINITY);

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
				rep[idx0]==-1 && rep[idx1]!=-1)
			x[rep[idx1]] = a(x[rep[idx1]], x[idx0]);
	}

	// fill-in the computed optima
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ij = j  * w + i;
		if (rep[ij]!=-1)
			x[ij] = x[rep[ij]];
	}

	//cleanup
	free(rep);
}


#ifndef OMIT_BDINT_MAIN
#define USE_BDINT_MAIN
#endif

#ifdef USE_BDINT_MAIN
#include <stdio.h>
#include "iio.h"
#include "smapa.h"
SMART_PARAMETER_SILENT(BDMAX,0)
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s in.tiff out.tiff\n", *v);
		//                          0 1       2
		return 1;
	}
	char *filename_in   = v[1];
	char *filename_out  = v[2];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	bdint(x, w, h, BDMAX()>0?fmaxf:fminf);

	iio_save_image_float(filename_out, x, w, h);

	return 0;
}
#endif//USE_BDINT_MAIN
