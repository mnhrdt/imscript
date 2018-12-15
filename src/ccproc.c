// This file contains the implemantation of the "ccproc" function
// It is a general tool to deal with connected components in images.
// Also, it contains the implementation of the "bfollow" function, to track
// the boundaries of these connected components in anticlockwise order.

#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include "abstract_dsf.c"
#include "xmalloc.c"

// equivalence relations of floats
//
// A function of this signature is the most important imput parameter.  It is
// used to decide whether neighboring pixels belong to the same connected
// component or not.
typedef int (*float_equivalence_relation_t)(float,float);


// three example equivalence relations are provided:

// wheter the two floats are equal
int float_equality(float a, float b)
{
	return a == b;
}

// whether two floats are equal, or are both null
int floatnan_equality(float a, float b)
{
	if (isnan(a) && isnan(b)) return true;
	return a == b;
}

// whether they are both nan or both non-nan
int float_eq_isnan(float a, float b)
{
	return isnan(a) == isnan(b);
}


// function to sort int vectors with increasing first component
static int compare_pair_increasing(const void *aa, const void *bb)
{
	const int *a = (const int *)aa;
	const int *b = (const int *)bb;
	return (*a > *b) - (*a < *b);
}

// function to sort int vectors with decreasing first component
static int compare_pair_decreasing(const void *aa, const void *bb)
{
	return -compare_pair_increasing(aa, bb);
}

// given an image of indexes, whether a given pixel is at a boundary
// of a region of pixels of the same index
static bool idx_boundaryingP(int *idxs, int w, int h, int idx)
{
	int i = idx % w;
	int j = idx / w;
	assert(j*w+i == idx);
	int my_idx = idxs[idx];
	if (i-1 >= 0 && my_idx != idxs[(j+0)*w + i-1]) return true;
	if (i+1 <  w && my_idx != idxs[(j+0)*w + i+1]) return true;
	if (j-1 >= 0 && my_idx != idxs[(j-1)*w + i+0]) return true;
	if (j+1 <  h && my_idx != idxs[(j+1)*w + i+0]) return true;
	return false;
}

static void swapi(int *t, int i, int j)
{
	int tmp = t[i];
	t[i] = t[j];
	t[j] = tmp;
}

// whether a point is inside the image domain
static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

// directions
//	0: ( 1, 0)
//	1: ( 0,-1)
//	2: (-1, 0)
//	3: ( 0, 1)

// get the neighbor of pixel (i,j) in the direction dir
static void get_neighbour(int *ni, int *nj, int i, int j, int dir)
{
	*ni = i;
	*nj = j;
	switch (dir) {
	case 0:	*ni += 1; break;
	case 1:	*nj -= 1; break;
	case 2:	*ni -= 1; break;
	case 3:	*nj += 1; break;
	}
}

// whether an edge (i,j,dir) is at the boundary of an eq-region
static bool pix_boundaryingP(float *x, int w, int h,
		float_equivalence_relation_t eq, int i, int j, int dir)
{
	assert(insideP(w, h, i, j));
	int ni, nj;
	get_neighbour(&ni, &nj, i, j, dir);
	if (!insideP(w, h, ni, nj))
		return true;
	float vij = x[ j*w+ i];
	float vnn = x[nj*w+ni];
	return !eq(vij, vnn);
}

// follow the closed boundary of a connected component containing point (i,j)
// note: the output array must be pre-allocated
// out_bd is a list of triplets (i,j,dir) conforming the boundary
int bfollow(int *out_bd, float *x, int w, int h,
		float_equivalence_relation_t eq,
		int i, int j)
{
	assert(i >= 0 && j >= 0 && i < w && j < h);

	// find a pixel that is outside, it will be used to begin the loop
	int ifirst = i;
	int jfirst = j;
	int dfirst = 1;
	while (!pix_boundaryingP(x, w, h, eq, ifirst, jfirst, dfirst))
		ifirst -= 1;

	// loop until the first boundary element is found again
	int out_n = 0;
	int dir = dfirst;
	i = ifirst;
	j = jfirst;
	do {
		assert(pix_boundaryingP(x, w, h, eq, i, j, dir));
		out_bd[3*out_n + 0] = i;
		out_bd[3*out_n + 1] = j;
		out_bd[3*out_n + 2] = dir;
		out_n += 1;
		switch(dir) { // counterclockwise track (could be a table also)
		case 0: if (pix_boundaryingP(x, w, h, eq, i, j, 1))
				dir = 1;
			else if (pix_boundaryingP(x, w, h, eq, i, j-1, 0))
				j -= 1;
			else if (pix_boundaryingP(x, w, h, eq, i+1, j-1, 3)) {
				i += 1; j -= 1; dir = 3;
			} else fail("bullshit 0\n");
			break;
		case 1: if (pix_boundaryingP(x, w, h, eq, i, j, 2))
				dir = 2;
			else if (pix_boundaryingP(x, w, h, eq, i-1, j, 1))
				i -= 1;
			else if (pix_boundaryingP(x, w, h, eq, i-1, j-1, 0)) {
				i -= 1; j -= 1; dir = 0;
			} else fail("bullshit 1\n");
			break;
		case 2: if (pix_boundaryingP(x, w, h, eq, i, j, 3))
				dir = 3;
			else if (pix_boundaryingP(x, w, h, eq, i, j+1, 2))
				j += 1;
			else if (pix_boundaryingP(x, w, h, eq, i-1, j+1, 1)) {
				i -= 1; j += 1; dir = 1;
			} else fail("bullshit 2\n");
			break;
		case 3: if (pix_boundaryingP(x, w, h, eq, i, j, 0))
				dir = 0;
			else if (pix_boundaryingP(x, w, h, eq, i+1, j, 3))
				i += 1;
			else if (pix_boundaryingP(x, w, h, eq, i+1, j+1, 2)) {
				i += 1; j += 1; dir = 2;
			} else fail("bullshit 3\n");
			break;
		}
	} while (i != ifirst || j != jfirst || dir != dfirst);
	return out_n;
}

// Compute the number of connected components of equivalent pixels in an image.
// If any of the "out_*" parameters is non-null, it is filled-in.
// The connected components are ordered by size.
// Note: connexity means 4-connexity.
//
// return value: N = the number of connected components
// out_size[i] = area of ith cc (0<=i<N)
// out_bdsize[i] = perimenter of ith cc (0<=i<N)
// out_all[p] = permutation of {0..w*h-1} where the cc are contiguous
// out_first[i] = first point of ith region on table out_all (0<=i<N)
// out_idx[p] = index of the region containing pixel p
//
int ccproc_old(
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
		eq = floatnan_equality;

	// create table of representatives
	int *rep = xmalloc(w * h * sizeof*rep);
	for (int i = 0; i < w*h; i++)
		rep[i] = i;
	adsf_begin(rep, w*h);

	// join equivalent neighbors (neighbors == 4-neighbors ALWAYS)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (i+1 < w && eq(x[p0], x[p1])) adsf_union(rep, w*h, p0, p1);
		if (j+1 < h && eq(x[p0], x[p2])) adsf_union(rep, w*h, p0, p2);
	}

	// canonicalize dsf (after this, the DSF is not changed anymore)
	for (int i = 0; i < w*h; i++)
		rep[i] = adsf_find(rep, w*h, i);
	//img_debug_int(rep, w, h, 1, "ccproc_rep.tiff");

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
	qsort(pair, r, 2*sizeof*pair, compare_pair_decreasing);
	for (int i = 0; i < r; i++)
	{
		int ir = pair[2*i + 1];
		tmpi[ir] = i;
		assert(ir == rep[ir]);
	}

	// fill-in table of sizes
	for (int i = 0; i < r; i++)
		out_size[i] = pair[2*i + 0];

	// fill-in image of representatives
	for (int i = 0; i < w*h; i++)
		out_idx[i] = tmpi[rep[i]];

	// reverse indexes
	int *ipair = xmalloc(2 * w * h * sizeof*ipair);
	for (int i = 0; i < w*h; i++)
	{
		ipair[2*i + 0] = tmpi[rep[i]];
		ipair[2*i + 1] = i;
	}
	qsort(ipair, w*h, 2*sizeof*ipair, compare_pair_increasing);
	for (int i = 0; i < w*h; i++)
		out_all[i] = ipair[2*i+1];

	// table of first indexes
	out_first[0] = 0;
	rcx = 1;
	for (int i = 1; i < w*h; i++)
		if (rep[out_all[i-1]] != rep[out_all[i]])
			out_first[rcx++] = i;
	assert(rcx == r);

	// verify consistency
	for (int i = 0; i < r; i++)
	for (int j = 0; j < out_size[i]; j++)
		assert(out_idx[out_all[out_first[i] + j]] == i);

	// reorder out_all so that boundarying points come first
	for (int i = 0; i < r; i++)
	{
		int *t = out_all + out_first[i];
		int a = 0;
		int b = out_size[i];
		while (a < b)
			if (idx_boundaryingP(out_idx, w, h, t[a]))
				a = a + 1;
			else
				swapi(t, a, --b);
		out_bdsize[i] = a;
	}

	// cleanup and exit
	free(ipair);
	free(pair);
	free(rep);
	free(tmpi);
	return r;
}

int ccproc_miniold(
		int *out_size,          // total size of each CC
		int *out_bdsize,        // boundary size of each CC
		int *out_all,           // array of all the indexes, connected
		int *out_first,         // array of first indexes of each CC
		int *out_idx,           // image with the indices of each region
		float *x, int w, int h, // input image
		float_equivalence_relation_t eq
	  )
{
	if (!eq) eq = floatnan_equality;          // set equivalence relation
	int *rep   = xmalloc(  w*h*sizeof*rep);   // representatives
	int *pair  = xmalloc(2*w*h*sizeof*pair);  // pairs
	int *tmpi  = xmalloc(  w*h*sizeof*tmpi);  // indexes
	int *ipair = xmalloc(2*w*h*sizeof*ipair); // inverse pairs
	for (int i = 0; i < w*h; i++) rep[i] = i; // set representatives
	adsf_begin(rep, w*h);                     // init dsf
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {  // for each pixel (i,j)
		int p0 = j*w + i;      // index of       (i,j)
		int p1 = j*w + i+1;    // index of pixel (i+1,j)
		int p2  = (j+1)*w + i; // index of pixel (i,j+1)
		if (i+1 < w && eq(x[p0], x[p1]))adsf_union(rep, w*h, p0, p1);
		if (j+1 < h && eq(x[p0], x[p2])) adsf_union(rep, w*h, p0, p2);
	}
	for (int i = 0; i < w*h; i++) rep[i] = adsf_find(rep, w*h, i);
	int R = 0; // number of regions (output value of this function)
	for (int i = 0; i < w*h; i++) if (rep[i] == i) R += 1;
	for (int i = 0; i < R; i++) pair[2*i + 0] = 0;
	int rcx = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i) {
			tmpi[i] = rcx;
			pair[2*rcx + 1] = i; // representative
			rcx += 1;
		}
	for (int i = 0; i < w*h; i++) pair[2*tmpi[rep[i]] + 0] += 1;
	qsort(pair, R, 2*sizeof*pair, compare_pair_decreasing);
	for (int i = 0; i <  R ; i++) tmpi[pair[2*i+1]] = i;
	for (int i = 0; i <  R ; i++) out_size[i] = pair[2*i + 0];
	for (int i = 0; i < w*h; i++) out_idx[i] = tmpi[rep[i]];
	for (int i = 0; i < w*h; i++) ipair[2*i + 0] = tmpi[rep[i]];
	for (int i = 0; i < w*h; i++) ipair[2*i + 1] = i;
	qsort(ipair, w*h, 2*sizeof*ipair, compare_pair_increasing);
	for (int i = 0; i < w*h; i++) out_all[i] = ipair[2*i+1];
	out_first[0] = 0;
	rcx = 1;
	for (int i = 1; i < w*h; i++)
		if (rep[out_all[i-1]] != rep[out_all[i]])
			out_first[rcx++] = i;
	for (int i = 0; i < R; i++) {
		int *t = out_all + out_first[i];
		int a = 0;
		int b = out_size[i];
		while (a < b)
			if (idx_boundaryingP(out_idx, w, h, t[a]))
				a = a + 1;
			else
				swapi(t, a, --b);
		out_bdsize[i] = a;
	}
	free(ipair);
	free(pair);
	free(rep);
	free(tmpi);
	return R;
}

static int compute_representatives(int *rep, float *x, int w, int h,
		float_equivalence_relation_t eq)
{
	// rep[i] = i
	adsf_begin(rep, w*h);

	// join equivalent neighbors (neighbors == 4-neighbors ALWAYS)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (i+1 < w && eq(x[p0], x[p1])) adsf_union(rep, w*h, p0, p1);
		if (j+1 < h && eq(x[p0], x[p2])) adsf_union(rep, w*h, p0, p2);
	}

	// canonicalize dsf (after this, the DSF is not changed anymore)
	for (int i = 0; i < w*h; i++)
		rep[i] = adsf_find(rep, w*h, i);

	// count connected components
	int r = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i)
			r += 1;
	return r;
}

static
void compute_pos_and_idx(int *pos, int *idx, int *rep, int w, int h, int r)
{
	for (int i = 0; i < r; i++)
		pos[2*i + 0] = 0; // initialize areas of each region

	int rcx = 0;
	for (int i = 0; i < w*h; i++)
		if (rep[i] == i)
		{
			idx[i] = rcx;
			pos[2*rcx + 1] = i; // representative
			rcx += 1;
		}
	assert(rcx == r);
	for (int i = 0; i < w*h; i++)
	{
		int j = rep[i];
		int t = idx[j];
		pos[2*t + 0] += 1;
	}
	qsort(pos, r, 2*sizeof*pos, compare_pair_decreasing);
	for (int i = 0; i < r; i++)
	{
		int ir = pos[2*i + 1];
		idx[ir] = i;
		assert(ir == rep[ir]);
	}
}

static void reverse_indexes(int *out_all, int *idx, int *rep, int n)
{
	int *ipair = xmalloc(2 * n * sizeof*ipair);

	// reverse indices
	for (int i = 0; i < n; i++)
	{
		ipair[2*i + 0] = idx[rep[i]];
		ipair[2*i + 1] = i;
	}
	qsort(ipair, n, 2*sizeof*ipair, compare_pair_increasing);

	// fill-in table of all indices
	for (int i = 0; i < n; i++)
		out_all[i] = ipair[2*i+1];

	free(ipair);
}

static void reorder_boundaries(int *bdsize,
	int *idx, int *all, int *first, int *size, int w, int h, int r)
{
	for (int i = 0; i < r; i++)
	{
		int *t = all + first[i];
		int a = 0;
		int b = size[i];
		// vector: t[0] t[1] ... t[b-1]
		while (a < b)
			if (idx_boundaryingP(idx, w, h, t[a]))
				a = a + 1;
			else
				swapi(t, a, --b);
		bdsize[i] = a;
	}
}

static
void assert_consistency(int *idx, int *all, int *size, int *first, int r)
{
	for (int i = 0; i < r; i++)
	for (int j = 0; j < size[i]; j++)
		assert(idx[all[first[i] + j]] == i);
}

// core function
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
		eq = floatnan_equality;

	// create table of representatives
	int *rep = xmalloc(w * h * sizeof*rep);
	int r = compute_representatives(rep, x, w, h, eq);

	// set the index/position correspondence
	int *pair = xmalloc(2 * r * sizeof*pair);
	int *tmpi = xmalloc(w*h*sizeof*tmpi);
	compute_pos_and_idx(pair, tmpi, rep, w, h, r);

	// reverse indexes and fill table with all
	reverse_indexes(out_all, tmpi, rep, w*h);

	// fill-in table of sizes
	for (int i = 0; i < r; i++)
		out_size[i] = pair[2*i + 0];

	// fill-in image of representatives
	for (int i = 0; i < w*h; i++)
		out_idx[i] = tmpi[rep[i]];

	// fill-in table of first indexes
	out_first[0] = 0;
	int cx = 1;
	for (int i = 1; i < w*h; i++)
		if (rep[out_all[i-1]] != rep[out_all[i]])
			out_first[cx++] = i;
	assert(cx == r);

	// verify consistency
	assert_consistency(out_idx, out_all, out_size, out_first, r);

	// reorder out_all so that boundarying points come first
	reorder_boundaries(out_bdsize, out_idx, out_all, out_first, out_size,
	                   w, h, r);

	// cleanup and exit
	free(pair);
	free(rep);
	free(tmpi);
	return r;
}

int ccproc_mini( // minified version of the function above
	int *o_size, int *o_bdsize, int *o_all, int *o_first, int *o_idx,
	float *x, int w, int h, float_equivalence_relation_t eq)
{
	if (!eq)
		eq = floatnan_equality;
	int rep[w*h];
	int r = compute_representatives(rep, x, w, h, eq);
	int pos[2*r];
	int idx[w*h];
	compute_pos_and_idx(pos, idx, rep, w, h, r);
	reverse_indexes(o_all, idx, rep, w*h);
	for (int i = 0; i < r; i++)
		o_size[i] = pos[2*i + 0];
	for (int i = 0; i < w*h; i++)
		o_idx[i] = idx[rep[i]];
	o_first[0] = 0;
	int cx = 1;
	for (int i = 1; i < w*h; i++)
		if (rep[o_all[i-1]] != rep[o_all[i]])
			o_first[cx++] = i;
	reorder_boundaries(o_bdsize, o_idx, o_all, o_first, o_size, w, h, r);
	return r;
}


#ifdef MAIN_BDFOLLOW
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 4)
		return fprintf(stderr, "usage:\n\t%s img i j\n", *v);
		//                                 0 1   2 3
	char *filename_in = v[1];
	int in_i = atoi(v[2]);
	int in_j = atoi(v[3]);

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	int *bd = xmalloc(3*4*w*h*sizeof*bd);

	int r = bfollow(bd, x, w, h, floatnan_equality, in_i, in_j);
	printf("r = %d\n", r);
	for (int i = 0; i < r; i++)
		printf("\t%d %d %d\n", bd[3*i+0], bd[3*i+1], bd[3*i+2]);

	return 0;
}
#endif//MAIN_BDFOLLOW

//#define MAIN_CCPROC
#ifdef MAIN_CCPROC
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
	int *out_size2 = xmalloc(w*h*sizeof*out_size);

	int r = ccproc(out_size, out_bdsize, out_all, out_first, out_idx,
			x, w, h,
			//float_eq_isnan
			floatnan_equality
			);

	fprintf(stderr, "ccproc returned %d\n", r);
	for (int i = 0; i < r; i++)
	{
		fprintf(stderr, "out_size[%d/%d] = %d (%d)\n",
				i, r, out_size[i], out_bdsize[i]);
		if (i > 100) break;
	}

	// verify consistency
	int totsize = 0;
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < out_size[i]; j++)
			x[out_all[out_first[i] + j]] = i;
		totsize += out_size[i];
	}
	assert(totsize == w*h);

	iio_write_image_int("ccproc_idx.tiff", out_idx, w, h);
	iio_write_image_int("ccproc_all.tiff", out_all, w, h);
	iio_write_image_float("ccproc_xxx.tiff", x, w, h);

	for (int i = 0; i < w*h; i++)
		x[i] = -1;
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < out_bdsize[i]; j++)
			x[out_all[out_first[i] + j]] = i;
	}
	iio_write_image_float("ccproc_yyy.tiff", x, w, h);


	for (int i = 0; i < w*h; i++)
		out_size2[i] = out_size[ out_idx[i] ] < 10;
	iio_write_image_int("ccproc_siz.tiff", out_size2, w, h);


	free(out_size);
	free(out_size2);
	free(out_bdsize);
	free(out_all);
	free(out_first);
	free(out_idx);
	free(x);
	return 0;
}
#endif

