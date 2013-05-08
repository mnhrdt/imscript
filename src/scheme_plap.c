#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// the type of a "getpixel" function
typedef double (*getpixel_operator)(double*,int,int,int,int);

// extrapolate by 0
inline static double getpixel_0(double *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w];
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
inline static double getpixel_1(double *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}


// global variable containing p
static double global_variable_containing_p = NAN;

// the type of a compact scheme
typedef double (*scheme9_t)(double [3][3]);

// u00 u10 u20
// u01 u11 u21
// u02 u12 u22

// the five-point scheme for the laplacian
static double scheme9_lap5(double u[3][3])
{
	return -4 * u[1][1] + u[1][0] + u[0][1] + u[1][2] + u[2][1];
}

// the compact nine-point scheme for the laplacian
static double scheme9_lap9(double u[3][3])
{
	return -6 * u[1][1] + u[1][0] + u[0][1] + u[1][2] + u[2][1]
		+ (u[0][0] + u[2][0] + u[0][2] + u[2][2])/2;
}


// a naif, symmetric, finite difference scheme for the infinity laplacian
static double scheme9_linf_naif(double u[3][3])
{
	// u00 u10 u20
	// u01 u11 u21
	// u02 u12 u22

	double uX = (u[2][1] - u[0][1])/2;
	double uY = (u[1][2] - u[1][0])/2;
	double uXX = -2*u[1][1] + u[2][1] + u[0][1];
	double uYY = -2*u[1][1] + u[1][2] + u[1][0];
	double uXY = -2*u[1][1] + u[1][0] + u[2][1] + u[1][2] + u[0][1]
		- u[0][2] - u[2][0];

	double r = uXX*uX*uX + uYY*uY*uY + uXY*uX*uY;
	return r;
}

// upwind five-point scheme for the infinity laplacian
static double scheme9_linf_mono4(double u[3][3])
{
	double min = INFINITY, max = -INFINITY;
	min = fmin(min, u[1][0]); max = fmax(max, u[1][0]);
	min = fmin(min, u[0][1]); max = fmax(max, u[0][1]);
	min = fmin(min, u[2][1]); max = fmax(max, u[2][1]);
	min = fmin(min, u[1][2]); max = fmax(max, u[1][2]);
	double r = min + max - 2*u[1][1];
	return r;
}

// a naif scheme for the p-laplacian
static double scheme9_lapp_naif(double u[3][3])
{
	// u00 u10 u20
	// u01 u11 u21
	// u02 u12 u22

	double l2 = scheme9_lap5(u);
	double linf = scheme9_linf_naif(u);

	double uX = (u[2][1] - u[0][1])/2;
	double uY = (u[1][2] - u[1][0])/2;
	double ngrad = hypot(uX, uY);

	double p = global_variable_containing_p;

	double r = pow(ngrad, p-2)*l2 + (p-2)*pow(ngrad, p-4)*linf;
	return r;
}

// an asymmetric divergence-form scheme for the p-laplacian
static double scheme9_lapp_div(double u[3][3])
{
	// g = grad(u)
	// (a,b) = |g|^(p-2) * g
	// r = div(a,b)

	// 00 10 20
	// 01 11 21
	// 02 12 22

	double p = global_variable_containing_p;

	double g11[2] = {u[2][1] - u[1][1], u[1][2] - u[1][1]};
	double g10[2] = {u[2][0] - u[1][0], u[1][1] - u[1][0]};
	double g01[2] = {u[1][1] - u[0][1], u[0][2] - u[0][1]};

	double npg11 = pow(hypot(g11[0], g11[1]), p);
	double npg10 = pow(hypot(g10[0], g10[1]), p);
	double npg01 = pow(hypot(g01[0], g01[1]), p);
	double a11 = npg11 * g11[0];
	double b11 = npg11 * g11[1];
	double a01 = npg01 * g01[0];
	double b10 = npg10 * g10[1];
	double r = a11 - a01 + b11 - b10;
	return r;
}

// the type of a "eval scheme" function
typedef double (*eval_scheme9_t)(double *,int,int,int,int,scheme9_t);

// evaluate a scheme at a point in an image, as is
static double eval_scheme9(double *x, int w, int h, int i, int j, scheme9_t s)
{
	getpixel_operator p = getpixel_1;

	double u[3][3];
	for (int ii = 0; ii < 3; ii++)
	for (int jj = 0; jj < 3; jj++)
		u[ii][jj] = p(x, w, h, i + ii - 1, j + jj - 1);

	double r = s(u);

	return r;
}

// evaluate a scheme at a point in an image, with symmetrization
static double eval_scheme9_sym(double *x, int w, int h, int i, int j, scheme9_t s)
{
	getpixel_operator p = getpixel_1;

	double u[4][3][3] = {
		{ {0,0,0}, {0,0,0}, {0,0,0} },
		{ {0,0,0}, {0,0,0}, {0,0,0} },
		{ {0,0,0}, {0,0,0}, {0,0,0} },
		{ {0,0,0}, {0,0,0}, {0,0,0} }
	};
	for (int ii = 0; ii < 3; ii++)
	for (int jj = 0; jj < 3; jj++)
	{
		double y = p(x, w, h, i + ii - 1, j + jj - 1);
		u[0][   ii ][   jj ] = y;
		u[1][ 2-ii ][   jj ] = y;
		u[2][   ii ][ 2-jj ] = y;
		u[3][ 2-ii ][ 2-jj ] = y;
	}

	double r = 0;
	for (int k = 0; k < 4; k++)
		r += s(u[k])/4;

	return r;
}

// fill an image using the provided scheme
static void fill_with_scheme(double *y, double *x, int w, int h,
		scheme9_t s, eval_scheme9_t es)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[j*w+i] = es(x, w, h, i, j, s);

}

#include "fail.c"

static void apply_plap(double *y, double *x, int w, int h, double p, char *s)
{
	global_variable_containing_p = p;

	if (false) ;
	else if (p == 2 && 0 == strcmp(s, "lap5"))
		fill_with_scheme(y,x,w,h, scheme9_lap5, eval_scheme9);
	else if (p == 2 && 0 == strcmp(s, "lap9"))
		fill_with_scheme(y,x,w,h, scheme9_lap9, eval_scheme9);
	else if (isinf(p) && 0 == strcmp(s, "mono4"))
		fill_with_scheme(y,x,w,h, scheme9_linf_mono4, eval_scheme9);
	else if (isinf(p) && 0 == strcmp(s, "naif"))
		fill_with_scheme(y,x,w,h, scheme9_linf_naif, eval_scheme9);
	else if (isinf(p) && 0 == strcmp(s, "naifsym"))
		fill_with_scheme(y,x,w,h, scheme9_linf_naif, eval_scheme9_sym);
	else if (isfinite(p) && 0 == strcmp(s, "naif"))
		fill_with_scheme(y,x,w,h, scheme9_lapp_naif, eval_scheme9);
	else if (isfinite(p) && 0 == strcmp(s, "naifsym"))
		fill_with_scheme(y,x,w,h, scheme9_lapp_naif, eval_scheme9_sym);
	else if (isfinite(p) && 0 == strcmp(s, "div"))
		fill_with_scheme(y,x,w,h, scheme9_lapp_div, eval_scheme9);
	else if (isfinite(p) && 0 == strcmp(s, "divsym"))
		fill_with_scheme(y,x,w,h, scheme9_lapp_div, eval_scheme9_sym);
	else
		fail("can not use scheme \"%s\" with p=%g", s, p);
}

#include "iio.h"
#include "xmalloc.c"

int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s {lap5|naif|mono4|div|fv} p [in [out]]\n", *v);
		//        0  1                       2  3   4
		return 1;
	}
	char *scheme_id = v[1];
	double p = atof(v[2]);
	char *infile = c > 3 ? v[3] : "-";
	char *outfile = c > 4 ? v[4] : "-";

	int w, h;
	double *x = iio_read_image_double(infile, &w, &h);
	double *y = xmalloc(w*h*sizeof*y);
	apply_plap(y, x, w, h, p, scheme_id);
	iio_save_image_double(outfile, y, w, h);
	free(x);
	return 0;
}
