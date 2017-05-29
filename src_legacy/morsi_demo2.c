// the simplest implementation of morphological operators

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "iio.h"


// utility function that returns a valid pointer to memory
static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new) {
		fprintf(stderr, "ERROR: out of memory when requesting "
			       "%zu bytes\n", size);
		abort();
	}
	return new;
}

// the type of the "getpixel" operators
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by nan
static float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return NAN;
	return x[i + j*w];
}


// a structuring element is a list E of integers
// E[0] = number of pixels
// E[1] = 0 (flags, not used yet)
// (E[2], E[3]) = position of the center
// (E[4], E[5]) = first pixel
// (E[6], E[7]) = second pixel
// ...

void morsi_erosion(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmin(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

void morsi_dilation(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = -INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmax(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

static float median(float *a, int n)
{
	if (n < 1) return NAN;
	if (n == 1) return *a;
	if (n == 2) return (a[0] + a[1])/2;
	qsort(a, n, sizeof*a, compare_floats);
	if (0 == n%2)
		return (a[n/2]+a[1+n/2])/2;
	else
		return a[n/2];
}

void morsi_median(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a[e[0]];
		int cx = 0;
		for (int k = 0; k < e[0]; k++)
		{
			float v = p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]);
			if (isfinite(v))
				a[cx++] = v;
		}
		y[j*w+i] = median(a, cx);
	}
}

//void morsi_opening(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_erosion(t, x, w, h, e);
//	morsi_dilation(y, t, w, h, e);
//	free(t);
//}
//
//void morsi_closing(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_dilation(t, x, w, h, e);
//	morsi_erosion(y, t, w, h, e);
//	free(t);
//}
//
//void morsi_gradient(float *y, float *x, int w, int h, int *e)
//{
//	float *a = xmalloc(w*h*sizeof*a);
//	float *b = xmalloc(w*h*sizeof*b);
//	morsi_erosion(a, x, w, h, e);
//	morsi_dilation(b, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = b[i] - a[i];
//	free(a);
//	free(b);
//}
//
//void morsi_igradient(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_erosion(t, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = x[i] - t[i];
//	free(t);
//}
//
//void morsi_egradient(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_dilation(t, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = t[i] - x[i];
//	free(t);
//}
//
//void morsi_laplacian(float *y, float *x, int w, int h, int *e)
//{
//	float *a = xmalloc(w*h*sizeof*a);
//	float *b = xmalloc(w*h*sizeof*b);
//	morsi_erosion(a, x, w, h, e);
//	morsi_dilation(b, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = (a[i] + b[i] - 2*x[i])/2;
//	free(a);
//	free(b);
//}
//
//void morsi_enhance(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_laplacian(t, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = x[i] - t[i];
//	free(t);
//}
//
//void morsi_strange(float *y, float *x, int w, int h, int *e)
//{
//	float *a = xmalloc(w*h*sizeof*a);
//	float *b = xmalloc(w*h*sizeof*b);
//	morsi_opening(a, x, w, h, e);
//	morsi_closing(b, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = b[i] - a[i];
//	free(a);
//	free(b);
//}
//
//void morsi_tophat(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_opening(t, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = x[i] - t[i];
//	free(t);
//}
//
//void morsi_bothat(float *y, float *x, int w, int h, int *e)
//{
//	float *t = xmalloc(w*h*sizeof*t);
//	morsi_closing(t, x, w, h, e);
//	for (int i = 0; i < w*h; i++)
//		y[i] = t[i] - x[i];
//	free(t);
//}

void morsi_all(
	float *o_ero, float *o_dil, float *o_med, float *o_ope, float *o_clo,
	float *o_grad, float *o_igrad, float *o_egrad,
	float *o_lap, float *o_enh, float *o_str,
	float *o_top, float *o_bot, float *x, int w, int h, int *e)
{
	int n = w*h, s = n*sizeof(float), own_ope=0, own_clo=0, own_lap=0, i;
	float *min = xmalloc(s);
	float *max = xmalloc(s);
	if ((o_top || o_str) && !o_ope) { o_ope = xmalloc(s); own_ope = 1; }
	if ((o_bot || o_str) && !o_clo) { o_clo = xmalloc(s); own_clo = 1; }
	if (o_enh && ! o_lap)           { o_lap = xmalloc(s); own_lap = 1; }

	morsi_erosion (min, x, w, h, e);
	morsi_dilation(max, x, w, h, e);
	if (o_ero)   for(i=0;i<n;i++) o_ero[i] = min[i];
	if (o_dil)   for(i=0;i<n;i++) o_dil[i] = max[i];
	if (o_med)   morsi_median(o_med, x, w, h, e);
	if (o_ope)   morsi_dilation(  o_ope,     min, w, h, e);
	if (o_clo)   morsi_erosion (  o_clo,     max, w, h, e);
	if (o_grad)  for(i=0;i<n;i++) o_grad[i]  = max[i] - min[i];
	if (o_igrad) for(i=0;i<n;i++) o_igrad[i] =   x[i] - min[i];
	if (o_egrad) for(i=0;i<n;i++) o_egrad[i] = max[i] -   x[i];
	if (o_top)   for(i=0;i<n;i++) o_top[i]   =     x[i] - o_ope[i];
	if (o_bot)   for(i=0;i<n;i++) o_bot[i]   = o_clo[i] -   x[i];
	if (o_str)   for(i=0;i<n;i++) o_str[i]   = o_clo[i] - o_ope[i];
	if (o_lap)   for(i=0;i<n;i++) o_lap[i]   = (max[i] + min[i] - 2*x[i])/2;
	if (o_enh)   for(i=0;i<n;i++) o_enh[i]   =     x[i] - o_lap[i];

	free(min); free(max);
	if (own_ope) free(o_ope);
	if (own_clo) free(o_clo);
	if (own_lap) free(o_lap);
}

static float *flinc(float *x, int nf)
{
	return x ? x + nf : x;
}

void morsi_all_split(
	float *o_ero, float *o_dil, float *o_med, float *o_ope, float *o_clo,
	float *o_grad, float *o_igrad, float *o_egrad,
	float *o_lap, float *o_enh, float *o_str,
	float *o_top, float *o_bot, float *x, int w, int h, int pd, int *e)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < pd; i++)
	{
		int n = i*w*h;
		morsi_all(flinc(o_ero,n), flinc(o_dil,n), flinc(o_med,n),
			flinc(o_ope,n), flinc(o_clo,n),
			flinc(o_grad,n), flinc(o_igrad,n), flinc(o_egrad,n),
			flinc(o_lap,n), flinc(o_enh,n), flinc(o_str,n),
			flinc(o_top,n), flinc(o_bot,n), flinc(x,n), w, h, e);
	}
}


static int *build_disk(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4, elen = 2*side*side+4;
	int *e = xmalloc(elen*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
	for (int j = -radius-1; j <= radius+1; j++)
		if (hypot(i,j) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = j;
			cx += 1;
		}
	assert(cx < side*side);
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_dysk(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4, elen = 2*side*side+4;
	int *e = xmalloc(elen*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
	for (int j = -radius-1; j <= radius+1; j++)
		if (hypot(i,j) < radius && hypot(i,j) >= radius-1) {
			e[2*cx+4] = i;
			e[2*cx+5] = j;
			cx += 1;
		}
	assert(cx < side*side);
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_hrec(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = 0;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_vrec(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = 0;
			e[2*cx+5] = i;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}


static int *build_structuring_element_from_string(char *name)
{
	static int cross[] = {5,0,  0,0, -1,0, 0,0, 1,0, 0,-1, 0,1 };
	static int square[] = {9,0, 0,0,
			-1,-1,-1,0,-1,1, 0,-1,0,0,0,1, 1,-1,1,0,1,1};
	int *r = NULL;
	if (0 == strcmp(name, "cross" )) r = cross;
	if (0 == strcmp(name, "square")) r = square;
	if (4 == strspn(name, "disk"  )) r = build_disk(atof(name+4));
	if (4 == strspn(name, "dysk"  )) r = build_dysk(atof(name+4));
	if (4 == strspn(name, "hrec"  )) r = build_hrec(atof(name+4));
	if (4 == strspn(name, "vrec"  )) r = build_vrec(atof(name+4));
	return r;
}

int main(int c, char **v)
{
	// process input arguments
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s element image\n", *v);
		//                          0 1       2
	int *structuring_element = build_structuring_element_from_string(v[1]);
	if (!structuring_element)
		return fprintf(stderr, "elements = cross, square, disk5 ...\n");
	char *filename_in = v[2];

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// alloc space for output images
	int n = w*h*pd;
	float *y = xmalloc(13*n*sizeof*y);

	// compute filters
	morsi_all_split(y+0*n, y+1*n, 0/*y+2*n*/, y+3*n, y+4*n, y+5*n, y+6*n,
			y+7*n, y+8*n, y+9*n, y+10*n, y+11*n, y+12*n,
			x, w, h, pd, structuring_element);

	// normalize the laplacian
	for(int i = 0; i < w*h*pd; i++)
		y[8*n+i] += 127;

	// save images
	struct { char *n; float *x; } t[13] = {
		{"001_erosion.png",   y+0*n },
		{"002_dilation.png",  y+1*n },
		{"003_median.png",    y+2*n },
		{"004_opening.png",   y+3*n },
		{"005_closing.png",   y+4*n },
		{"006_gradient.png",  y+5*n },
		{"007_igradient.png", y+6*n },
		{"008_egradient.png", y+7*n },
		{"009_laplacian.png", y+8*n },
		{"010_enhance.png",   y+9*n },
		{"011_oscill.png",    y+10*n},
		{"012_tophat.png",    y+11*n},
		{"013_bothat.png",    y+12*n}
	};

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < 13; i++)
		if (t[i].x)
			iio_write_image_float_split(t[i].n, t[i].x, w, h, pd);

	// cleanup and exit
	free(x); free(y);
	return 0;
}
