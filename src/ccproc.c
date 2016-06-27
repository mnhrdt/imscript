#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include "abstract_dsf.c"
#include "xmalloc.c"

typedef int (*float_equivalence_relation_t)(float,float);

int float_equality(float a, float b)
{
	return a == b;
}

int floatnan_equality(float a, float b)
{
	if (isnan(a) && isnan(b)) return true;
	return a == b;
//	uint32_t *aa = (void*)&a;
//	uint32_t *bb = (void*)&b;
//	return aa == bb;
}

static int compare_pair(const void *aa, const void *bb)
{
	const int *a = (const int *)aa;
	const int *b = (const int *)bb;
	return (*a > *b) - (*a < *b);
}

static int ncompare_pair(const void *aa, const void *bb)
{
	const int *a = (const int *)aa;
	const int *b = (const int *)bb;
	return (*a < *b) - (*a > *b);
}

static bool boundaryingP(int *idxs, int w, int h, int idx)
{
	int i = idx % w;
	int j = idx / w;
	assert(j*w+i == idx);
	int my_idx = idxs[idx];
	if (i-1 >= 0 && my_idx != idxs[(j+0)*w + i-1]) return false;
	if (i+1 <  w && my_idx != idxs[(j+0)*w + i+1]) return false;
	if (j-1 >= 0 && my_idx != idxs[(j-1)*w + i+0]) return false;
	if (j+1 <  h && my_idx != idxs[(j+1)*w + i+0]) return false;
	return true;
}

static void swapi(int *t, int i, int j)
{
	int tmp = t[i];
	t[i] = t[j];
	t[j] = t[i];
}

// compute the number of connected components of equivalent pixels in an image
// if any of the "out_*" parameters is non-null, it is filled-in
// the connected components are ordered by size
//
// return value: N = the number of connected components
// out_size[i] = area of ith cc (0<=i<N)
// out_bdsize[i] = perimenter of ith cc (0<=i<N)
// out_all[p] = permutation of {0..w*h-1} where the cc are contiguous
// out_first[i] = first point of ith region on table out_all (0<=i<N)
// out_idx[p] = index of the region containing pixel p
//
int ccproc(
		int *out_size,          // total size of each CC
		int *out_bdsize,        // boundary size of each CC
		int *out_all,           // array of all the indexes, connected
		int *out_first,         // array of first indexes of each CC
		int *out_idx,           // image with the indices of each region
		float *x, int w, int h, // input image
		float_equivalence_relation_t eq
	  )
{
	// pre-setup
	if (!eq)
		eq = float_equality;

	// create table of representatives
	int *rep = xmalloc(w * h * sizeof*rep);
	for (int i = 0; i < w*h; i++)
		rep[i] = i;
	adsf_begin(rep, w*h);

	// join equivalent neighbors (neighbors == 4-neighbors ALWAYS)
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (eq(x[p0], x[p1])) adsf_union(rep, w*h, p0, p1);
		if (eq(x[p0], x[p2])) adsf_union(rep, w*h, p0, p2);
	}

	// canonicalize dsf (after this, the DSF is not changed anymore)
	for (int i = 0; i < w*h; i++)
		rep[i] = adsf_find(rep, w*h, i);

	// count connected components
	int r = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i)
			r += 1;

	// for now, assume that everywhing is there
	assert(out_size);
	assert(out_bdsize);
	assert(out_all);
	assert(out_first);
	assert(out_idx);

	// set the index/position correspondence
	int *pair = xmalloc(2 * r * sizeof*pair);
	int *tmpi = xmalloc(w*h*sizeof*tmpi);
	for (int i = 0; i < r; i++)
		pair[2*i + 0] = 0; // initialize areas of each region
	int rcx = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i)
		{
			tmpi[i] = rcx;
			pair[2*rcx + 1] = i; // representative
			rcx += 1;
		}
	assert(rcx == r);
	for (int i = 0; i < w*h; i++)
	{
		int j = rep[i];
		int t = tmpi[j];
		pair[2*t + 0] += 1;
	}
	qsort(pair, r, 2*sizeof*pair, ncompare_pair);
	for (int i = 0; i < r; i++)
	{
		int ir = pair[2*i + 1];
		tmpi[ir] = i;
		assert(ir == rep[ir]);
	}

	// fill-in table of sizes
	if (out_size)
		for (int i = 0; i < r; i++)
			out_size[i] = pair[2*i + 0];


	// fill-in image of representatives
	if (out_idx)
		for (int i = 0; i < w*h; i++)
			out_idx[i] = tmpi[rep[i]];

	// reverse indexes
	int *ipair = xmalloc(2 * w * h * sizeof*ipair);
	for (int i = 0; i < w*h; i++)
	{
		ipair[2*i + 0] = tmpi[rep[i]];
		ipair[2*i + 1] = i;
	}
	qsort(ipair, w*h, 2*sizeof*ipair, compare_pair);
	for (int i = 0; i < w*h; i++)
		out_all[i] = ipair[2*i+1];

	// table of first indexes
	out_first[0] = 0;
	rcx = 1;
	for (int i = 1; i < w*h; i++)
		if (rep[out_all[i-1]] != rep[out_all[i]])
			out_first[rcx++] = i;
	fprintf(stderr, "rcx = %d\n", rcx);
	fprintf(stderr, "r = %d\n", r);
	assert(rcx == r);

	// verify consistence
	for (int i = 0; i < r; i++)
	for (int j = 0; j < out_size[i]; j++)
		assert(out_idx[out_all[out_first[i] + j]] == i);

	// reorder out_all so that boundarying points come first
	for (int i = 0; i < r; i++)
	{
		int *ti = out_all + out_first[i];
		int topo = 0;
		for (int j = 0; j < out_size[i]; j++)
		{
			int idx = ti[j];
			if (!boundaryingP(out_idx, w, h, idx))
			{
				swapi(ti, topo, j);
				topo += 1;
			}
		}
		out_bdsize[i] = topo;
	}

	// cleanup and exit
	free(pair);
	free(rep);
	free(tmpi);
	return r;
}

#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s img.png\n", *v);
	char *filename_in = v[1];

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	int *out_size = xmalloc(w*h*sizeof*out_size);
	int *out_bdsize = xmalloc(w*h*sizeof*out_size);
	int *out_all = xmalloc(w*h*sizeof*out_size);
	int *out_first = xmalloc(w*h*sizeof*out_size);
	int *out_idx = xmalloc(w*h*sizeof*out_size);

	int r = ccproc(out_size, out_bdsize, out_all, out_first, out_idx,
			x, w, h, floatnan_equality);

	fprintf(stderr, "ccproc returned %d\n", r);
	for (int i = 0; i < r; i++)
	{
		fprintf(stderr, "out_size[%d/%d] = %d (%d)\n",
				i, r, out_size[i], out_bdsize[i]);
		if (i > 100) break;
	}

	// verify consistence
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < out_size[i]; j++)
			x[out_all[out_first[i] + j]] = i;
	}

	iio_save_image_int("ccproc_idx.tiff", out_idx, w, h);
	iio_save_image_int("ccproc_all.tiff", out_all, w, h);
	iio_save_image_float("ccproc_xxx.tiff", x, w, h);

	for (int i = 0; i < w*h; i++)
		x[i] = -1;
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < out_bdsize[i]; j++)
			x[out_all[out_first[i] + j]] = i;
	}
	iio_save_image_float("ccproc_yyy.tiff", x, w, h);


	return 0;
}
