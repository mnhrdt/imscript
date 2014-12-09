// simulator of affine satellite acquisitions
// INPUT:
// 	1. height map
// 	2. texture map
// 	3. projection matrix
// 	4. desired output size
//
// OUTPUT:
// 	1. rendered texture

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fail.c"

static void invert_projection(double iA[8], double A[8])
{
	// 
	//  / i \     / a  b \     / l \     / r \         / p \
	// |     | = |        | * |     | + |     | * h + |     |
	//  \ j /     \ c  d /     \ t /     \ s /         \ q /
	//
	double a = A[0], b = A[1], r = A[2], p = A[3];
	double c = A[4], d = A[5], s = A[6], q = A[7];

	double det = a * d - b * c;
	double ia =  d / det;
	double ib = -b / det;
	double ic = -c / det;
	double id =  a / det;
	double ir = -(ia * r + ib * s);
	double is = -(ic * r + id * s);
	double ip = -(ia * p + ib * q);
	double iq = -(ic * p + id * q);

	iA[0] = ia; iA[1] = ib; iA[2] = ir; iA[3] = ip;
	iA[4] = ic; iA[5] = id; iA[6] = is; iA[7] = iq;
}

static void apply_projection(double y[3], double A[8], double x[3])
{
	y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2] + A[3];
	y[1] = A[4] * x[0] + A[5] * x[1] + A[6] * x[2] + A[7];
	y[2] = x[2];
}

#include "getpixel.c"

static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

float bicubic_interpolation(float *img, int w, int h, float x, float y)
{
	x -= 1;
	y -= 1;

	getsample_operator p = getsample_0;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, 1, ix + i, iy + j, 0);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}

void bicubic_interpolation_vec(float *result,
		float *img, int w, int h, int pd, float x, float y)
{
	x -= 1;
	y -= 1;

	getsample_operator p = getsample_0;

	int ix = floor(x);
	int iy = floor(y);
	for (int l = 0; l < pd; l++) {
		float c[4][4];
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
		float r = bicubic_interpolation_cell(c, x - ix, y - iy);
		result[l] = r;
	}
}


typedef double (objective_function_t)(double,void*);

static double bisection_1d(objective_function_t *f, double a, double b, void *e)
{
	double fa = f(a, e);
	double fb = f(b, e);
	if (fa * fb >= 0)
		return NAN;

	double c = ( a + b ) / 2;
	double fc = f(c, e);

	while (fabs(a - b) > 1e-6) {
		if (fc * fa >= 0) {
			a = c;
			fa = fc;
		} else {
			b = c;
			fb = fc;
		}
		c = ( a + b ) / 2;
		fc = f(c, e);
	}

	return c;
}


struct bisection_state {
	float *heights;
	int w, h;
	double *L;
	double ij[3];
};

double objective_height_side(double h, void *ee)
{
	// 0. The height "h" parametrizes positions on the line of view.
	struct bisection_state *e = ee;
	double pij[3] = {e->ij[0], e->ij[1], h};

	// 1. Given a height "h", compute a 3d point "p"
	double p[3];
	apply_projection(p, e->L, pij);

	// 2. Project this point into the plane "h=0"
	assert(p[2] == h);
	p[2] = 0;

	// 3. Compute the elevation at this position
	float elev = bicubic_interpolation(e->heights, e->w, e->h, p[0], p[1]);

	// 4. Substract this elevation from h
	return h - elev;
}

static void raytrace(double out[3],
		double L[8], float *heights, int w, int h, double ij[2])
{
	double base[3] = {ij[0], ij[1], 0};
	apply_projection(base, L, ij);
	double bh = bicubic_interpolation(heights, w, h, base[0], base[1]);

	struct bisection_state e[1];
	e->w = w;
	e->h = h;
	e->heights = heights;
	e->L = L;
	e->ij[0] = ij[0];
	e->ij[1] = ij[1];

	double hh = bisection_1d(objective_height_side, 0, 2*bh+1, e);
	double intersection[3] = {ij[0], ij[1], hh};
	apply_projection(out, L, intersection);


	if (ij[0] == 250 && ij[1] == 250) {
		fprintf(stderr, "raytrace check\n");
		fprintf(stderr, "\tij = %g %g\n", ij[0], ij[1]);
		fprintf(stderr, "\tbase = %g %g\n", base[0], base[1]);
		fprintf(stderr, "\tintersection = %g %g %g\n", intersection[0], intersection[1], intersection[2]);
		fprintf(stderr, "\tout = %g %g %g\n", out[0], out[1], out[2]);

	//	double h0 = objective_height_side(0, e);
	//	double h1 = objective_height_side(out[2], e);
	//	fprintf(stderr, "\th0,h1 = %g %g\n", h0, h1);
	//	//double hh = bisection_1d(objective_height_side, 0, out[2], e);
	//	double vhh = objective_height_side(hh, e);
	//	fprintf(stderr, "\thh, vhh = %g %g\n", hh, vhh);
	//	out[2] = hh;
	}
}

void satproj(float *out, int ow, int oh,
		double P[8], float *heights, float *colors,
		int w, int h, int pd)
{
	// P = projection
	// L = localisation
	//
	double L[8];
	invert_projection(L, P);

	fprintf(stderr, "satproj %d %d => %d %d (%d)\n", w, h, ow, oh, pd);
	fprintf(stderr, "projection P = %g %g %g %g  %g %g %g %g\n",
			P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
	fprintf(stderr, "projection L = %g %g %g %g  %g %g %g %g\n",
			L[0], L[1], L[2], L[3], L[4], L[5], L[6], L[7]);

	// fill each point of the output image with the appropriate color
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		double ij[3] = {i, j, 0}, rs[3] = {i, j, 0};
		raytrace(rs, L, heights, w, h, ij);
		float *to = out + pd * (ow * j + i);
		bicubic_interpolation_vec(to, colors, w, h, pd, rs[0], rs[1]);
	}
}



//double fun(double x, void *e)
//{
//	double r = x*x-2;
//	fprintf(stderr, "f(%.30f)=%g\n", x, r);
//	return r;
//}
//int main()
//{
//	double x = bisection_1d(fun, 0, 2, NULL);
//	printf("%lf\n", x);
//	return 0;
//}

static void center_projection(double P[8], double cx, double cy)
{
	// point "c" must be fixed at h = 0
	P[3] = cx - P[0]*cx - P[1]*cy;
	P[7] = cy - P[4]*cx - P[5]*cy;
}

#define MAIN_SATPROJ
#ifdef MAIN_SATPROJ
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
#include "parsenumbers.c"


static void read_n_doubles_from_string(double *out, char *string, int n)
{
	for (int i = 0; i < n; i++)
		out[i] = 0;

	int no;
	double *buf = NULL;
	FILE *f = fopen(string, "r");
	if (f) {
		buf = read_ascii_doubles(f, &no);
		fclose(f);
	} else {
		buf = alloc_parse_doubles(n, string, &no);
	}

	if (no > n) no = n;
	for (int i = 0; i < no; i++)
		out[i] = buf[i];
	free(buf);
}

#include <stdbool.h>
#include "pickopt.c"

int main(int c, char *v[])
{
	// input arguments
	bool do_center = pick_option(&c, &v, "c", NULL);
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s heights.tiff colors.png P.txt ow oh out.png\n", *v);
		//        0 1            2          3     4  5  6
		return 1;
	}
	char *fname_heights = v[1];
	char *fname_colors  = v[2];
	char *fname_pmatrix = v[3];
	int out_w = atoi(v[4]);
	int out_h = atoi(v[5]);
	char *fname_output  = v[6];

	// read input images and matrices
	int w[2], h[2], pd;
	float *heights = iio_read_image_float(fname_heights, w, h);
	float *colors  = iio_read_image_float_vec(fname_colors, w+1, h+1, &pd);
	double P[8] = {0};
	read_n_doubles_from_string(P, fname_pmatrix, 8);
	if (w[0] != w[1] || h[0] != h[1])
		fail("color and heights size mismatch");

	fprintf(stderr, "heights %d %d\n", w[0], h[0]);
	fprintf(stderr, "colors %d %d %d \n", w[1], h[1], pd);

	// allocate space for output
	float *out = xmalloc(out_w * out_h * pd * sizeof*out);

	// perform centering, if necessary
	if (do_center)
		center_projection(P, *w/2, *h/2);

	// run simulator
	satproj(out, out_w, out_h, P, heights, colors, *w, *h, pd);

	// save output
	iio_save_image_float_vec(fname_output, out, out_w, out_h, pd);

	// cleanup and exit
	free(heights);
	free(colors);
	free(out);
	return 0;
}
#endif//MAIN_SATPROJ
