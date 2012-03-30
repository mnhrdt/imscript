#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
	//extension_operator_float p = extend_float_image_by_zero;
	extension_operator_float p = extend_float_image_constant;
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
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ubar[j*w+i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}

static float sqr(float x)
{
	return x * x;
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
	//float maxdiff = 0;
	long double l2diff = 0;
	for (int i = 0; i < w*h; i++) {
		float t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		float newu = ubar[i] - Ex[i] * t;
		float newv = vbar[i] - Ey[i] * t;
		//if (fabs(newu - u[i]) > maxdiff) maxdiff = fabs(newu-u[i]);
		//if (fabs(newv - v[i]) > maxdiff) maxdiff = fabs(newv-v[i]);
		l2diff += sqr(newu-u[i])+sqr(newv-v[i]);
		u[i] = newu;
		v[i] = newv;
	}
	free(ubar);
	free(vbar);
	//return maxdiff < epsilon;
	return sqrt(l2diff/(w*h)) < epsilon;
}

void hs(float *u, float *v, float *a, float *b, int w, int h,
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

int hs_stopping(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha, float eps)
{
	//fprintf(stderr, "HSS N=%d a=%g e=%g\n", niter, alpha, eps);
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
	//fprintf(stderr, "HSS ran %d\n", i);
	free(gx);
	free(gy);
       	free(gt);
	return i;
}

#ifndef OMIT_MAIN
#include "iio.h"
int main(int argc, char *argv[])
{
	if (argc != 6 && argc != 7)
		exit(fprintf(stderr, "usage:\n\t%s niter alpha a b f\n",*argv));
	int niter = atoi(argv[1]);
	float alpha = atof(argv[2]);
	float epsilon = argc == 7 ? atof(argv[3]) : NAN;
	char *filename_a = argv[argc-3];
	char *filename_b = argv[argc-2];
	char *filename_f = argv[argc-1];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		exit(fprintf(stderr, "input images size mismatch\n"));
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	if (isfinite(epsilon)) {
		int nit = hs_stopping(u, v, a, b, w, h, niter, alpha, epsilon);
		fprintf(stderr, "ran %d iterations\n", nit);
	} else
		hs(u, v, a, b, w, h, niter, alpha);
	float *f = xmalloc(w * h * 2 * sizeof(float));
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_save_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
