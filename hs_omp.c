#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "iio.h"

static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "out of memory\n"));
	return new;
}

typedef float (*extension_operator_float)(float*, int, int, int, int);

static float extend_float_image_by_zero(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i > w-1 || j > h-1)
		return 0;
	else
		return x[j*w+i];
}

static float extend_float_image_constant(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static void compute_input_derivatives(float *Ex, float *Ey, float *Et,
		float *a, float *b, int w, int h)
{
	extension_operator_float p = extend_float_image_by_zero;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		//Ey[j*w+i] = (1.0/4) * ( p(a,w,h, i, j+1) - p(a,w,h, i,j-1)
		//		+p(b,w,h, i, j+1) - p(b,w,h, i,j-1));
		//Ex[j*w+i] = (1.0/4) * ( p(a,w,h, i+1, j) - p(a,w,h, i-1,j)
		//		+p(b,w,h, i+1, j) - p(b,w,h, i-1,j));
		//Et[j*w+i] = ( p(b,w,h, i, j) - p(a,w,h, i,j));
		Ey[j*w+i] = (1.0/4) * ( p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j));
		Ex[j*w+i] = (1.0/4) * ( p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1));
		Et[j*w+i] = (1.0/4) * ( p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1));
	}
}

static void compute_bar(float *ubar, float *u, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;
#pragma omp parallel for
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ubar[j*w+i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}

static void hs_iteration(float *u, float *v,
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

static bool hs_iteration_stopping(float *u, float *v,
		float *Ex, float *Ey, float *Et, int w, int h, float alpha,
		float epsilon)
{
	float *ubar = xmalloc(w * h * sizeof(float));
	float *vbar = xmalloc(w * h * sizeof(float));
	compute_bar(ubar, u, w, h);
	compute_bar(vbar, v, w, h);
	double scale = (19.0*INT_MAX)/(20+w*h);
	int mm = 0;
#pragma omp parallel for reduction(|:mm)
	for (int i = 0; i < w*h; i++) {
		float t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		float newu = ubar[i] - Ex[i] * t;
		float newv = vbar[i] - Ey[i] * t;
		float du = fabs(newu - u[i]);
		float dv = fabs(newv - v[i]);
		u[i] = newu;
		v[i] = newv;
		int mmnew = scale * (du > dv ? du : dv);
		//fprintf(stderr, "mmnew %x\n", mmnew);
		mm |= mmnew;
	}
	float mms = mm / scale;
	//fprintf(stderr, "scale =  %g\n", scale);
	//fprintf(stderr, "mm %o\n", mm);
	//fprintf(stderr, "mms %g\n", mms);
	free(ubar);
	free(vbar);
	return mms < epsilon;
}

static void hs(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (int i = 0; i < w*h; i++)
		u[i] = v[i] = 0;
	for (int i = 0; i < niter; i++)
		hs_iteration(u, v, gx, gy, gt, w, h, alpha);
	free(gx);
	free(gy);
       	free(gt);
}

static int hs_stopping(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha, float eps)
{
	int i;
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (i = 0; i < w*h; i++)
		u[i] = v[i] = 0;
	for (i = 0; i < niter; i++)
		if (hs_iteration_stopping(u, v, gx, gy, gt, w, h, alpha, eps))
			break;
	free(gx);
	free(gy);
       	free(gt);
	return i;
}

int main(int argc, char *argv[])
{
	if (argc != 8)
		exit(fprintf(stderr, "usage:\n\t"
			"%s nprocs niter alpha epsilon a b f\n",*argv));
	     //           0 1      2      3     4      5 6 7
	int nprocs = atoi(argv[1]);
	int niter = atoi(argv[2]);
	float alpha = atof(argv[3]);
	float epsilon = atof(argv[4]);
	char *filename_a = argv[5];
	char *filename_b = argv[6];
	char *filename_f = argv[7];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		exit(fprintf(stderr, "input images size mismatch\n"));
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	if(nprocs > 0) omp_set_num_threads(nprocs);
	int nit = hs_stopping(u, v, a, b, w, h, niter, alpha, epsilon);
	fprintf(stderr, "ran %d iterations\n", nit);
	float *f = xmalloc(w * h * 2 * sizeof(float));
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_save_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
