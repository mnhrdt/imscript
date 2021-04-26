#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#include "xmalloc.c"

struct floatintpos { float f; int i; int p[2]; } ;
static int compare_struct_floatintpos(const void *aa, const void *bb)
{
	const struct floatintpos *a = (const struct floatintpos *)aa;
	const struct floatintpos *b = (const struct floatintpos *)bb;
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
                                                 int i, int j, int ii, int jj,
						 float W)
{
	int nW = winpatch_length(pd, W);
	float px[nW], py[nW];
	getpatch(px, x, w, h, pd, i, j, W);
	getpatch(py, y, w, h, pd, ii, jj, W);
	return patchdistance(px, py, nW);
}

#define OMIT_BLUR_MAIN
#include "blur.c"

static void get_reference_image(float *, float **, int, int, int, int);

void silly_lucky_regions(float **y, float **x, int n, int w, int h, int pd,
		float W, float S, float R, float *init_img)
{
	// 1. compute or get the reference image
	// 2. compute the S-blurred version of each input frame
	// 3. for each pixel location
	//      3.1. compute the vector of W-patch-distances between
	//           the reference image and each frame
	//      3.2. sort the vector of distances, with their frame indices
	//      3.3. write the corresponding sorted values into the output video

	// 1. compute the average image
	float *mx = xmalloc(w*h*pd*sizeof*mx);
	if (init_img)
		for (int i = 0; i < w*h*pd; i++)
			mx[i] = init_img[i];
	else
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

	// 3.-1 (build a table of neighbors)
	int nn = 0, RR = R;
	if (RR < 0) RR = 0;
	if (RR > 30) RR = 30;
	int neigs[(2*RR+1)*(2*RR+1)][2];
	for (int i = -RR; i <= RR; i++)
	for (int j = -RR; j <= RR; j++)
	{
		neigs[nn][0] = i;
		neigs[nn][1] = j;
		nn += 1;
	}
	assert(nn == (2*RR+1)*(2*RR+1));

	// 3. for each pixel location
	fprintf(stderr, "processing lines...\n");
	time_t timeref = time(NULL);
	for (int j = 0; j < h; j++) {
		if (0==j%10) {
			int timtim = time(NULL) - timeref;
			fprintf(stderr, "\tline %d/%d {%d}\n", j, h-1, timtim);
			//fprintf(stderr, "\tline %d/%d\n", j,h-1);
		}
	for (int i = 0; i < w; i++)
	{
	//      3.1. compute the vector of winsize-patch-distances between
	//          the average image and each frame
		struct floatintpos v[n*nn];
		int cx = 0;
		for (int k = 0; k < n; k++)
		for (int l = 0; l < nn; l++)
		{
			int ii = i + neigs[l][0];
			int jj = j + neigs[l][1];
			if (ii >= 0 && jj >= 0 && ii < w && jj < h)
			{
				v[cx].i = k;
				v[cx].p[0] = ii;
				v[cx].p[1] = jj;
				v[cx].f = wpatch_distance(mx, bx[k], w, h, pd,
						i, j, ii, jj, W);
				cx += 1;
			}
		}
		assert(cx >= n);
		assert(cx < n*nn);

	//      3.2. sort the vector of distances, with their frame indices
		qsort(v, cx, sizeof*v, compare_struct_floatintpos);

	//      3.3. write the corresponding sorted values into the output video
		//int idx = j*w + i;
		for (int k = 0; k < n; k++)
			for (int l = 0; l < pd; l++)
			{
				int idx = j*w + i;
				int idx2 = v[k].p[1]*w + v[k].p[0];
				y[k][idx*pd+l] = x[v[k].i][idx2*pd+l];
			}
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
		fprintf(stderr, "computing the average image...\n");
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
	if (c != 8 && c != 9) {
		fprintf(stderr, "usage:\n\t"
			"%s wize ker rad inpat outpat first last [init]\n", *v);
		//        0 1    2   3   4     5      6     7     8
		return 1;
	}
	float winsize = atof(v[1]);
	float kersigma = atof(v[2]);
	float searchrad = atof(v[3]);
	char *filepattern_in = v[4];
	char *filepattern_out = v[5];
	int idx_first = atoi(v[6]);
	int idx_last = atoi(v[7]);
	char *filename_init = c > 8 ? v[8] : NULL;

	int n = idx_last - idx_first + 1;
	if (n < 3 || n > 1000) fail("weird n = %d\n", n);

	// read input images
	fprintf(stderr, "reading input images...\n");
	int w=0, h=0, pd = 0;
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

	// read initialization image, if supplied
	float *init_img = NULL;
	if (filename_init) {
		int ww, hh, ppdd;
		init_img = iio_read_image_float_vec(filename_init,
							&ww, &hh, &ppdd);
		if (w != ww || h != hh || pd != ppdd)
			fail("init image size mismatch");
	}

	// do stuff
	fprintf(stderr, "computing stuff...\n");
	silly_lucky_regions(y, x, n, w, h, pd,
			winsize, kersigma, searchrad,
			init_img);

	// save output images
	fprintf(stderr, "saving output images...\n");
	for (int i = 0; i < n; i++)
	{
		int idx = idx_first + i;
		char filename_out[FILENAME_MAX];
		snprintf(filename_out, FILENAME_MAX, filepattern_out, idx);
		iio_write_image_float_vec(filename_out, y[i], w, h, pd);
	}

	return 0;
}
