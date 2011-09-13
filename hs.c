// implementation of Horn & Schunck method for two frames

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"


// auxiliary function to abort the program with an error message
static void abort_with_message(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#ifdef NDEBUG
	exit(-1);
#else//NDEBUG
	exit(*(int *)0x43);
#endif//NDEBUG
}

// auxiliary function to request memory
// it always returns a valid block
static void *xmalloc(size_t size)
{
	if (size == 0)
		abort_with_message("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
	{
		double sm = size / (0x100000 * 1.0);
		abort_with_message("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)", size, sm);
	}
	return new;
}



typedef float (*extension_operator_float)(float*,int,int,int,int);

static float extend_float_image_by_zero(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0 || j < 0 || i > w-1 || j > h-1)
		return 0;
	else
		return x[j][i];
}

static float extend_float_image_constant(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j][i];
}





// Compute the spatial and temporal derivatives of the input data.
// This is copied verbatim from Section 7 of the original paper.
// The only difference is that the order of i and j is swapped.
static void compute_input_derivatives(float *gx, float *gy, float *gt,
		float *a, float *b, int w, int h)
{
	float (*Ex)[w] = (void*)gx;
	float (*Ey)[w] = (void*)gy;
	float (*Et)[w] = (void*)gt;

	// TODO: probably is better to extend images by constant, but the
	// original paper says zero, so we do it.
	extension_operator_float p = extend_float_image_by_zero;

	// TODO (?): re-write the following scheme by separating it into
	// 1) convolution with constant 2x2x2 kernel
	// 2) forward diffeences
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		Ey[j][i] = (1.0/4) * (
				p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j)
				);
		Ex[j][i] = (1.0/4) * (
				p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1)
				);
		Et[j][i] = (1.0/4) * (
				p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1)
				);
	}
}

// Compute the local averaging of a velocity field.
// This is copied verbatim from Section 8 of the original paper.
static void compute_bar(float *ubar, float *u, int w, int h)
{
	float (*ub)[w] = (void*)ubar;

	extension_operator_float p = extend_float_image_constant;

	// NOTE: this is the bottleneck of the code
	// the following loop amounts for 80% of the running time
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ub[j][i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));

	// ALTERNATIVE FOR OPTIMIZATION:
	//
	//  1     / 1 2 1 |
	// --- *  | 2 4 2 |   -   |3|
	// 12     \ 1 2 1 /
	//          ^^^^^
	//        separable

}

// Run one iteration of the optimization method.
// This is copied verbatim from Section 12 of the original paper.
static void hs_iteration(float *u, float *v, float *a, float *b,
		float *Ex, float *Ey, float *Et, int w, int h, float alpha)
{
	float *ubar = xmalloc(w * h * sizeof(float));
	float *vbar = xmalloc(w * h * sizeof(float));

	compute_bar(ubar, u, w, h);
	compute_bar(vbar, v, w, h);

	for (int i = 0; i < w*h; i++) {
		float t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		u[i] = ubar[i] - Ex[i] * t;
		v[i] = vbar[i] - Ey[i] * t;
	}

	free(ubar);
	free(vbar);
}

// Run "niter" iterations of the optimization method with parameter "alpha"
static void hs(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));

	compute_input_derivatives(gx, gy, gt, a, b, w, h);

	iio_save_image_float("/tmp/Ex", gx, w, h);
	iio_save_image_float("/tmp/Ey", gy, w, h);
	iio_save_image_float("/tmp/Et", gt, w, h);

	// initial solution
	for (int i = 0; i < w*h; i++)
		u[i] = v[i] = 0;

	// iterate
	for (int i = 0; i < niter; i++)
		hs_iteration(u, v, a, b, gx, gy, gt, w, h, alpha);

	free(gx);
	free(gy);
	free(gt);
}

// command-line interface
int main(int argc, char *argv[])
{
	if (argc != 6) {
		fprintf(stderr, "usage:\n\t%s niter alpha a b f\n", *argv);
		return EXIT_FAILURE;
	}
	int niter = atoi(argv[1]);
	float alpha = atoi(argv[2]);
	char *filename_a = argv[3];
	char *filename_b = argv[4];
	char *filename_f = argv[5];

	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		abort_with_message("two frames sizing mismatch");

	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	hs(u, v, a, b, w, h, niter, alpha);

	float *f = xmalloc(w * h * 2 * sizeof(float));
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_save_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
