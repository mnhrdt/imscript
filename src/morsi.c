#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "xmalloc.c"
//static void *xmalloc(size_t size)
//{
//	void *new = malloc(size);
//	if (!new) {
//		fprintf(stderr, "ERROR: out of memory when requesting "
//			       "%zu bytes\n", size);
//		abort();
//	}
//	return new;
//}

typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

// extrapolate by nan
static float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return NAN;
	return x[i + j*w];
}

// extrapolate by nearest value
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
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

void morsi_rank(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int cx = 0;
		float u = p(x,w,h, i, j);
		for (int k = 0; k < e[0]; k++)
		{
			float v = p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]);
			if (isfinite(v))
				cx += v < u;
		}
		y[j*w+i] = cx;
	}
}

void morsi_opening(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	morsi_dilation(y, t, w, h, e);
	free(t);
}

void morsi_closing(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_dilation(t, x, w, h, e);
	morsi_erosion(y, t, w, h, e);
	free(t);
}

void morsi_gradient(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_erosion(a, x, w, h, e);
	morsi_dilation(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = b[i] - a[i];
	free(a);
	free(b);
}

void morsi_igradient(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_egradient(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_dilation(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = t[i] - x[i];
	free(t);
}

void morsi_laplacian(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_erosion(a, x, w, h, e);
	morsi_dilation(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = (a[i] + b[i] - 2*x[i])/2;
	free(a);
	free(b);
}

void morsi_enhance(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_laplacian(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_blur(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_laplacian(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] + t[i];
	free(t);
}

void morsi_oscillation(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_opening(a, x, w, h, e);
	morsi_closing(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = b[i] - a[i];
	free(a);
	free(b);
}

void morsi_tophat(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_opening(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_bothat(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_closing(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = t[i] - x[i];
	free(t);
}

void morsi_iblur(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = (x[i] + t[i])/2;
	free(t);
}

void morsi_eblur(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_dilation(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = (x[i] + t[i])/2;
	free(t);
}

void morsi_cblur(float *y, float *x, int w, int h, int *e)
{
	float *t1 = xmalloc(w*h*sizeof*t1);
	float *t2 = xmalloc(w*h*sizeof*t1);
	morsi_erosion(t1, x, w, h, e);
	morsi_dilation(t2, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = 0.5*x[i] + 0.25*t1[i] + 0.25*t2[i];
	free(t1);
	free(t2);
}


void morsi_all(
	float *o_ero, float *o_dil, float *o_ope, float *o_clo,
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


int *build_disk(float radius)
{
	if (!(radius >1)) return NULL;
	//fprintf(stderr, "building a disk of radius %g\n", radius);
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
	//fprintf(stderr, "building a dysk of radius %g\n", radius);
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
	//fprintf(stderr, "building a hrec of radius %g\n", radius);
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
	//fprintf(stderr, "building a vrec of radius %g\n", radius);
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

static int *build_drec(float radius)
{
	if (!(radius >1)) return NULL;
	//fprintf(stderr, "building a drec of radius %g\n", radius);
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = i;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_Drec(float radius)
{
	if (!(radius >1)) return NULL;
	//fprintf(stderr, "building a drec of radius %g\n", radius);
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = -i;
			e[2*cx+5] = i;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

#ifndef OMIT_MORSI_MAIN
#define MORSI_TEST_MAIN
#endif

#ifdef MORSI_TEST_MAIN
#include <string.h>
#include "iio.h"
static char *help_string_name     = "morsi";
static char *help_string_version  = "morsi 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "simple gray-scale morphology";
static char *help_string_usage    = "usage:\n\t"
"morsi {cross|square|disk5...} {erosion|dilation...} [in [out]]";
static char *help_string_long     =
"Morsi applies a gray-scale morphological operator to an image.\n"
"\n"
"You have to specify a structuring element and a morphological operation.\n"
"There are no defaults.\n"
"\n"
"Usage: morsi ELEMENT OPERATION in out\n"
"   or: morsi ELEMENT OPERATION in > out\n"
"   or: cat in | morsi ELEMENT OPERATION > out\n"
"\n"
"Elements:\n"
" cross        a cross of 5 pixels\n"
" square       a square of 3x3 pixels\n"
" diskR        a discrete disk of radius R (e.g., disk3.1)\n"
" dyskR        the boundary of a discrete disk of radius R\n"
" hrecW        a rectangle of size Wx1\n"
" vrecH        a rectangle of size 1xH\n"
" drecX        a 45-degree diagonal of X-pixels\n"
" DrecX        a -45-degree diagonal of X-pixels\n"
"\n"
"Operations:\n"
" erosion      min of neighboring pixels (smooth and darken)\n"
" dilation     max of neighboring pixels (smooth and brighten)\n"
" median       median of neighboring pixels (smooth)\n"
" opening      erosion and dilation (remove small dark spots)\n"
" closing      dilation and erosion (remove small bright spots)\n"
" gradient     dilation minus erosion (find centered boundaries)\n"
" igradient    image minus erosion (inner boundaries)\n"
" egradient    dilation minus image (outer boundaries)\n"
" laplacian    difference between outer and inner grads (signed boundaries)\n"
" enhance      image minus laplacian (sharpen boundaries and details)\n"
" blur         image plus laplacian (centered smoothing, removes deatils)\n"
" tophat       image minus opening (select bright blobs)\n"
" bothat       closing minus image (select dark blobs)\n"
" oscillation  closing minus opening (measure local oscillation)\n"
" iblur        average of image and its erosion (blur towards dark)\n"
" eblur        average of image and its dilation (blur towards light)\n"
" cblur        average of iblur and eblur (detail-preserving smoothing)\n"
"\n"
"Examples:\n"
" morsi cross erosion i.png o.png    Erode by a \"cross\" structuring element\n"
" morsi disk4.2 median i.png o.png   Median filtering of radius 4.2\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above

int main_morsi(int c, char **v)
{
	// process "help" arguments
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// data for available structuring elements
	int cross[] = {5,0,  0,0, -1,0, 0,0, 1,0, 0,-1, 0,1 };
	int square[] = {9,0, 0,0, -1,-1,-1,0,-1,1, 0,-1,0,0,0,1, 1,-1,1,0,1,1};

	// process input arguments
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
				"%s element operation [in [out]]\n", *v);
		//                0 1       2          3   4
		return 1;
	}
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";
	int *structuring_element = NULL;
	if (0 == strcmp(v[1], "cross" )) structuring_element = cross;
	if (0 == strcmp(v[1], "square")) structuring_element = square;
	if (4==strspn(v[1],"disk"))structuring_element=build_disk(atof(v[1]+4));
	if (4==strspn(v[1],"dysk"))structuring_element=build_dysk(atof(v[1]+4));
	if (4==strspn(v[1],"hrec"))structuring_element=build_hrec(atof(v[1]+4));
	if (4==strspn(v[1],"vrec"))structuring_element=build_vrec(atof(v[1]+4));
	if (4==strspn(v[1],"drec"))structuring_element=build_drec(atof(v[1]+4));
	if (4==strspn(v[1],"Drec"))structuring_element=build_Drec(atof(v[1]+4));
	if (!structuring_element) {
		fprintf(stderr, "elements = cross, square ...\n");
		return 1;
	}
	void (*operation)(float*,float*,int,int,int*) = NULL;
	if (0 == strcmp(v[2], "erosion"    )) operation = morsi_erosion;
	if (0 == strcmp(v[2], "dilation"   )) operation = morsi_dilation;
	if (0 == strcmp(v[2], "median"     )) operation = morsi_median;
	if (0 == strcmp(v[2], "rank"       )) operation = morsi_rank;
	if (0 == strcmp(v[2], "opening"    )) operation = morsi_opening;
	if (0 == strcmp(v[2], "closing"    )) operation = morsi_closing;
	if (0 == strcmp(v[2], "gradient"   )) operation = morsi_gradient;
	if (0 == strcmp(v[2], "igradient"  )) operation = morsi_igradient;
	if (0 == strcmp(v[2], "egradient"  )) operation = morsi_egradient;
	if (0 == strcmp(v[2], "laplacian"  )) operation = morsi_laplacian;
	if (0 == strcmp(v[2], "enhance"    )) operation = morsi_enhance;
	if (0 == strcmp(v[2], "blur"       )) operation = morsi_blur;
	if (0 == strcmp(v[2], "oscillation")) operation = morsi_oscillation;
	if (0 == strcmp(v[2], "tophat"     )) operation = morsi_tophat;
	if (0 == strcmp(v[2], "bothat"     )) operation = morsi_bothat;
	if (0 == strcmp(v[2], "iblur"      )) operation = morsi_iblur;
	if (0 == strcmp(v[2], "eblur"      )) operation = morsi_eblur;
	if (0 == strcmp(v[2], "cblur"      )) operation = morsi_cblur;
	if (!operation) {
		fprintf(stderr, "operations = erosion, dilation, opening...\n");
		return 1;
	}

	// prepare input and output images
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = malloc(w*h*pd*sizeof*y);

	// compute
	if (pd == 1)
		operation(y, x, w, h, structuring_element);
	else
		for (int k = 0; k < pd; k++)
			operation(y+k*w*h, x+k*w*h, w, h, structuring_element);

	// save result
	iio_write_image_float_split(filename_out, y, w, h, pd);

	// cleanup
	free(x);
	free(y);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_morsi(c, v); }
#endif

#endif//MORSI_TEST_MAIN
