#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


//#include "xmalloc.c"

//// utility function that always returns a valid pointer to memory
//static void *xmalloc(size_t n)
//{
//	void *new = malloc(n);
//	if (!new)
//	{
//		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
//		exit(1);
//	}
//	return new;
//}
#include "conjugate_gradient.c"

// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[3]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[3] = xmalloc(w*h*3*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			mask[cx][2] = j*w+i;
			cx += 1;
		}
	assert(cx == nmask);

	*out_nmask = nmask;
	return mask;
}

static void invert_mask(int *im, int (*mask)[3], int nmask, int w, int h)
{
	for (int i = 0; i < w*h; i++)
		im[i] = -1;
	for (int i = 0; i < nmask; i++)
	{
		int idx = mask[i][2];// + w*mask[i][1];
		assert(idx >= 0);
		assert(idx < w*h);
		im[idx] = i;
	}
}



struct cglap_state {
	int w, h, (*mask)[3], nmask, *invmask;
	float *boundary_data;

};

typedef float (*fancy_getpixel_operator)(double*x,void*,int,int);

static float masked_getpixel_1(double *x, void *ee, int i, int j)
{
	struct cglap_state *e = ee;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= e->w) i = e->w-1;
	if (j >= e->h) j = e->h-1;
	int idx = j*e->w + i;
	assert(idx >= 0);
	assert(idx < e->w * e->h);
	//float *y = e->invmask[idx] < 0 ? e->original_image : x;
	//return y[idx];
	int pidx = e->invmask[idx];
	if (pidx < 0) {
		//assert(isfinite(e->original_image[idx]));
		return e->boundary_data[idx];
	} else {
		assert(pidx < e->nmask);
		//assert(isnan(e->original_image[idx]));
		return x[pidx];
	}
}

static double evaluate_laplacian_at(double *x, int p, void *ee)
{
	struct cglap_state *e = ee;
	assert(p >= 0);
	assert(p < e->nmask);
	int i = e->mask[p][0];
	int j = e->mask[p][1];

	fancy_getpixel_operator g = masked_getpixel_1;

	double r = -4 * g(x, ee, i  , j  )
		     + g(x, ee, i+1, j  )
		     + g(x, ee, i  , j+1)
		     + g(x, ee, i-1, j  )
		     + g(x, ee, i  , j-1);
	return r;
}

static void minus_laplacian_operator(double *y, double *x, int n, void *ee)
{
	for (int p = 0; p < n; p++)
		y[p] = -evaluate_laplacian_at(x, p, ee);
}

#include "smapa.h"
SMART_PARAMETER(CG_MAXIT,-1)
SMART_PARAMETER(CG_EPS,-1)

void harmonic_extension_steps(float *out, float *in, int w, int h)
{
	// build list of masked pixels
	int nmask, (*mask)[3] = build_mask(&nmask, in, w, h);
	int *invmask = xmalloc(w*h*sizeof(int));
	invert_mask(invmask, mask, nmask, w, h);

	// define the linear map A=laplacian_operator
	struct cglap_state e[1];
	e->w = w;
	e->h = h;
	e->mask = mask;
	e->invmask = invmask;
	e->nmask = nmask;
	e->boundary_data = xmalloc(w*h*sizeof(double));
	for (int i = 0; i < w*h; i++)
		e->boundary_data[i] = in[i];

	// fill-in the independent term b
	double *b = xmalloc(nmask * sizeof(double));
	double *tmp = xmalloc(nmask * sizeof(double));
	for (int p = 0; p < nmask; p++)
		tmp[p] = 0;
	minus_laplacian_operator(b, tmp, nmask, e);
	free(tmp);
	for (int p = 0; p < nmask; p++)
		b[p] *= -1;
	for (int i = 0; i < w*h; i++)
		if (e->invmask[i] < 0)
			e->boundary_data[i] = 0;

	// compute the solution
	double *solution = xmalloc(nmask * sizeof(double));
	double *initialization = xmalloc(nmask * sizeof(double));
	for (int i = 0; i < nmask; i++)
		initialization[i] = 0;
	//conjugate_gradient(solution, minus_laplacian_operator, b, nmask, e);
	int cg_maxit = CG_MAXIT() >= 0 ? CG_MAXIT() : nmask;
	float cg_eps = CG_EPS() >= 0 ? CG_EPS() : 1e-6;
	fancy_conjugate_gradient(solution, minus_laplacian_operator, b, nmask,
					e, initialization, cg_maxit, cg_eps);


	// copy the solution to its place
	for (int i = 0; i < w*h; i++)
		out[i] = in[i];
	for (int p = 0; p < nmask; p++) {
		if (nmask < 33)
			fprintf(stderr, "sol[%d] = %g\n", p, solution[p]);
		out[mask[p][2]] = solution[p];
	}

	free(e->boundary_data);
	free(solution);
	free(mask);
	free(invmask);
	free(b);
}

#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 4) {
		fprintf(stderr, "usage:\n\t"
			"%s in.png mask.png out.png\n", *argv);
		//        0 1      2        3
		return 1;
	}
	char *filename_in = argv[1];
	char *filename_mask = argv[2];
	char *filename_out = argv[3];

	int w[2], h[2], pd;
	float *in = iio_read_image_float(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	float *out = xmalloc(*w**h*pd*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			in[i] = NAN;

	harmonic_extension_steps(out, in, *w, *h);

	iio_save_image_float(filename_out, out, *w, *h);

	return 0;
}
