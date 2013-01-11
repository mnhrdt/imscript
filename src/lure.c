#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "xmalloc.c"

struct floatint { float f; int i; } ;
static int compare_struct_floatint(const void *aa, const void *bb)
{
	const struct floatint *a = (const struct floatint *)aa;
	const struct floatint *b = (const struct floatint *)bb;
	return (a->f - b->f) - (a->f < b->f);
}

// vector length of a patch of given witdth
static int winpatch_length(int pd, float W)
{
	int wi = W;
	int wis = 2*wi+1;
	int r = pd*wis*wis;
	static bool firsttime = true;
	if (firsttime) {
		firsttime = false;
		fprintf(stderr, "winpatch_length(%d, %g) = %d\n", pd, W, r);
	}
	return r;
}

#include "getpixel.c"

// fill a patch at a given location
static void getpatch(float *p, float *x, int w, int h, int pd,
                                         int i, int j, float W)
{
	int rad = W;
	int cx = 0;
	for (int di = -rad; di <= rad; di++)
	for (int dj = -rad; dj <= rad; dj++)
	for (int l = 0; l < pd; l++)
		p[cx++] = getsample_0(x, w, h, pd, i+di, j+dj, l);
		//p[cx++] = x[(w*j+i)*pd+l];
	assert(cx == winpatch_length(pd, W));
}

// distance between two patches
static float patchdistance(float *p, float *q, int n)
{
	// TODO (maybe) : use a fancier distance?
	float r = 0;
	for (int i = 0; i < n; i++)
		r = hypot(r, p[i] - q[i]);
	return r;
}

// compute the distance between two patches of size "W"
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

static void get_reference_image(float *, float **, int, int, int, int);

void silly_lucky_region(float **y, float **x, int n, int w, int h, int pd,
		float W, float S)
{
	// 1. compute the average image
	// 2. compute the S-blurred version of each input frame
	// 3. for each pixel location
	//      3.1. compute the vector of W-patch-distances between
	//           the average image and each frame
	//      3.2. sort the vector of distances, with their frame indices
	//      3.3. write the corresponding sorted values into the output video

	// 1. compute the average image
	fprintf(stderr, "computing the average image...\n");
	float *mx = xmalloc(w*h*pd*sizeof*mx);
	get_reference_image(mx, x, n, w, h, pd);

	// 2. compute the sigma-blurred version of each input frame
	fprintf(stderr, "computing the blurry input frames...\n");
	float *bx[n], blurparams[1] = {S};
	for (int i = 0; i < n; i++)
		bx[i] = xmalloc(w * h * pd * sizeof*bx[0]);
	for (int i = 0; i < n; i++)
		blur_2d(bx[i], x[i], w, h, pd, "gaussian", blurparams, 1);
	//	for (int j = 0; j < w*h*pd; j++)
	//		bx[i][j] = x[i][j];  // TODO: do the blur

	// 3. for each pixel location
	fprintf(stderr, "processing lines...\n");
	for (int j = 0; j < h; j++) {
		if (0==j%10)
			fprintf(stderr, "\tline %d/%d\n", j,h-1);
	for (int i = 0; i < w; i++)
	{
	//      3.1. compute the vector of winsize-patch-distances between
	//          the average image and each frame
		struct floatint v[n];
		for (int k = 0; k < n; k++)
		{
			v[k].i = k;
			v[k].f = wpatch_distance(mx, bx[k], w, h, pd, i, j, W);
		}

	//      3.2. sort the vector of distances, with their frame indices
		qsort(v, n, sizeof*v, compare_struct_floatint);

	//      3.3. write the corresponding sorted values into the output video
		int idx = j*w + i;
		for (int k = 0; k < n; k++)
			for (int l = 0; l < pd; l++)
				y[k][idx*pd+l] = x[v[k].i][idx*pd+l];
	}
	}

	// free temporary memory
	free(mx);
	for (int i = 0; i < n; i++)
		free(bx[i]);
}

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"

static void get_reference_image(float *mx,
		float **x, int n, int w, int h, int pd)
{
	char *fname_refimg = getenv("LURE_REFIMG");
	if (fname_refimg) { // use the provided file
		int rw, rh, rd;
		float *t = iio_read_image_float_vec(fname_refimg, &rw,&rh, &rd);
		if (!t)
			fail("could not open refimg \"%s\"", fname_refimg);
		if (rw != w || rh != h || rd != pd)
			fail("reference image size mismatch");
		for (int i = 0; i < w*h*pd; i++)
			mx[i] = t[i];
		free(t);
	} else { // compute the average
		for (int i = 0; i < w*h*pd; i++)
		{
			long double m = 0;
			for (int j = 0; j < n; j++)
				m += x[j][i];
			mx[i] = m/n;
		}
	}
}


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
	if (n < 3 || n > 1000) fail("weird n = %d\n", n);

	// read input images
	fprintf(stderr, "reading input images...\n");
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

	// allocate output images
	for (int i = 0; i < n; i++)
		y[i] = xmalloc(w * h * pd * sizeof*y[0]);

	// do stuff
	fprintf(stderr, "computing stuff...\n");
	silly_lucky_region(y, x, n, w, h, pd, winsize, kersigma);

	// save output images
	fprintf(stderr, "saving output images...\n");
	for (int i = 0; i < n; i++)
	{
		int idx = idx_first + i;
		char filename_out[FILENAME_MAX];
		snprintf(filename_out, FILENAME_MAX, filepattern_out, idx);
		iio_save_image_float_vec(filename_out, y[i], w, h, pd);
	}

	return 0;
}
