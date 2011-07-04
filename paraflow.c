// paraflow:
// Find a parametric motion model between two images by nonlinear
// multidimensional minimization


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_multimin.h>

#include "iio.h"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#include "fragments.c"



static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static float getsample(float *fx, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	float (*x)[w][pd] = (void*)fx;
	return x[j][i][l];
	//return x[(i+j*w)*pd + l];
}

static void bilinear_interpolation_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q)
{
	int ip = p;
	int iq = q;
	FORL(pd) {
		float a = getsample(x, w, h, pd, ip  , iq  , l);
		float b = getsample(x, w, h, pd, ip+1, iq  , l);
		float c = getsample(x, w, h, pd, ip  , iq+1, l);
		float d = getsample(x, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

//static void invflow(float *ou, float *flo, float *pin, int w, int h, int pd)
//{
//	float (*out)[w][pd] = (void*)ou;
//	float (*in)[w][pd] = (void*)pin;
//	float (*flow)[w][2] = (void*)flo;
//
//	FORJ(h) FORI(w) {
//		float p[2] = {i + flow[j][i][0], j + flow[j][i][1]};
//		float result[pd];
//		bilinear_interpolation_at(result, pin, w, h, pd, p[0], p[1]);
//		FORL(pd)
//			out[j][i][l] = result[l];
//	}
//}

static void affine_map(double y[2], double A[6], double x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

//static void invert_affinity(double invA[6], double A[6])
//{
//	double a, b, c, d, p, q;
//	a=A[0]; b=A[1]; p=A[2];
//	c=A[3]; d=A[4]; q=A[5];
//	double det = a*d - b*c;
//	invA[0] = d;
//	invA[1] = -b;
//	invA[2] = b*q-d*p;
//	invA[3] = -c;
//	invA[4] = a;
//	invA[5] = c*p-a*q;
//	FORI(6) invA[i] /= det;
//}

static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	//fprintf(stderr, "icell %g %g %g %g (%g %g)\n", a,b,c,d,x,y);
	switch(method) {
	//case 0: return interpolate_nearest(a, b, c, d, x, y);
	//case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	default: error("caca de vaca");
	}
	return -1;
}

static void general_interpolate(float *result,
		float *px, int w, int h, int pd, float p, float q,
		int m) // method
{
	float (*x)[w][pd] = (void*)px;
	if (p < 0 || q < 0 || p+1 >= w || q+1 >= h) {
		FORL(pd) result[l] = 0;
	} else {
		int ip = floor(p);
		int iq = floor(q);
		FORL(pd) {
			//float a = getsample(x, w, h, pd, ip  , iq  , l);
			//float b = getsample(x, w, h, pd, ip  , iq+1, l);
			//float c = getsample(x, w, h, pd, ip+1, iq  , l);
			//float d = getsample(x, w, h, pd, ip+1, iq+1, l);
			float a = x[iq][ip][l];
			float b = x[iq+1][ip][l];
			float c = x[iq][ip+1][l];
			float d = x[iq+1][ip+1][l];
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			//fprintf(stderr, "p%g q%g ip%d iq%d a%g b%g c%g d%g l%d v%g\n", p, q, ip, iq, a, b, c, d, l, v);
			result[l] = v;
		}
	}
}

static void apply_affinity(float *py, float *x, int w, int h, int pd,
		double A[6])
{
	float (*y)[w][pd] = (void*)py;
	//double invA[6]; invert_affinity(invA, A);
	FORJ(h) FORI(w) {
		double p[2] = {i, j}, q[2];
		affine_map(q, A, p);
		float val[pd];
		general_interpolate(val, x, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

//#include "vvector.h"
//static void invert_homography(double invH[9], double H[9])
//{
//	double h[3][3] = { {H[0], H[1], H[2]},
//			{H[3], H[4], H[5]},
//			{H[6], H[7], H[8]}};
//	double det;
//	double ih[3][3];
//	INVERT_3X3(ih, det, h);
//	FORI(9) invH[i] = ih[0][i];
//}

static void projective_map(double y[2], double H[9], double x[2])
{
	double z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
}

static void apply_homography(float *py, float *x, int w, int h, int pd,
			double H[9])
{
	float (*y)[w][pd] = (void*)py;
	//double invH[9]; invert_homography(invH, H);
	FORJ(h) FORI(w) {
		double p[2] = {i, j}, q[2];
		projective_map(q, H, p);
		float val[pd];
		general_interpolate(val, x, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}


//static double solvecubicspecial(double a, double b)
//{
//	long double x;
//	long double r;
//	if (a < 0) {
//		long double p = 1/a;
//		long double q = -b/a;
//		int k = 1;
//		long double cosarg = acos((3*q)/(2*p)*sqrt(-3/p))/3-k*2*M_PI/3;
//		r = 2*sqrt(-p/3) * cos(cosarg);
//		return r;
//	} else {
//		x = cbrt(sqrt((27*a*b*b+4)/a)/(2*sqrt(27)*a)+b/(2*a));
//		r = x-1/(3*a*x);
//		return r;
//	}
//}
//
//static double invertparabolicdistortion(double a, double xp)
//{
//	return solvecubicspecial(a, xp);
//}

// X' = X + a*X*X*X
static double parabolicdistortion(double a, double x)
{
	return x + a*x*x*x;
}

static void apply_radialpol(float *py, float *x, int w, int h, int pd,
		double *param)
{
	float (*y)[w][pd] = (void*)py;
	double c[2] = {param[0], param[1]};
	double a = param[2];
	fprintf(stderr, "radialpol (%g,%g)[%g]\n", c[0], c[1], a);
	FORJ(h) FORI(w) {
		double p[2] = {i, j}, q[2];
		double r = hypot(p[0]-c[0], p[1]-c[1]);
		double R = parabolicdistortion(a, r);
		//double ir = invertparabolicdistortion(a, r);
		// q = c + ir*(p-c)
		//q[0] = c[0] + (ir/r) * (p[0] - c[0]);
		//q[1] = c[1] + (ir/r) * (p[1] - c[1]);
		if (r > 0) {
			q[0] = c[0] + (R/r)*(i - c[0]);
			q[1] = c[1] + (R/r)*(j - c[1]);
		} else {
			q[0] = 0;
			q[1] = 0;
		}
		float val[pd];
		general_interpolate(val, x, w, h, pd, q[0], q[1], 2);
		FORL(pd) {
			y[j][i][l] = val[l];
		}
	}
}

static void apply_parametric_invflow(float *y, float *x, int w, int h, int pd,
		char *model_id, double *param, int nparam)
{
	FORI(w*h*pd) y[i] = -42;
	if (false) { ;
	} else if (0 == strcmp(model_id, "traslation")) {
		assert(nparam == 2);
		double a[6] = {1, 0, param[0], 0, 1, param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "euclidean")) {
		assert(nparam == 3);
		double a[6] = {cos(param[2]), sin(param[2]), param[0],
			-sin(param[2]), cos(param[2]), param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "similar")) {
		assert(nparam == 4);
		double r = param[3];
		double a[6] = {r*cos(param[2]), r*sin(param[2]), param[0],
			-r*sin(param[2]), r*cos(param[2]), param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "affine")) {
		assert(nparam == 6);
		double a[6]; FORI(6) a[i] = param[i];
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "projective")) {
		assert(nparam == 8);
		double a[9]; FORI(8) a[i] = param[i];
		a[8] = 1;
		apply_homography(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "cparabolic")) {
		assert(nparam == 1);
		double a[3] = {w/2.0, h/2.0, param[0]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "parabolic")) {
		assert(nparam == 3);
		double a[3] = {param[0], param[1], param[2]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "real1")) {
		error("not yet implemented");
	} else if (0 == strcmp(model_id, "real2")) {
		error("not yet implemented");
	} else error("unrecognized flow model %s", model_id);
}

//
// HERE BE DRAGONS (GSL multidimensional minimization)
//

struct problem_data {
	// image pair
	int w, h, pd;
	float *x, *y;

	// model metadata
	char *model_id;
};

static double objective_function_l2(const gsl_vector *v, void *pp)
{
	struct problem_data *p = pp;
	float *bp = xmalloc(p->w * p->h * p->pd * sizeof*bp);
	double point[v->size];
	FORI((int)v->size)
		point[i] = gsl_vector_get(v, i);
	apply_parametric_invflow(bp, p->x, p->w, p->h, p->pd, p->model_id,
			point, v->size);
	free(bp);
	return 0;
}


//
// END OF GSL STUFF
//

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

static void do_stuff(float *x, float *y, int w, int h, int pd,
		char *model_id, char *error_id,
		double *param, int nparams)
{
	float *iy = xmalloc(w*h*pd*sizeof*iy);
	apply_parametric_invflow(iy, y, w, h, pd, model_id, param, nparams);
	iio_save_image_float_vec("/tmp/merdota.tiff", iy, w, h, pd);
	free(iy);
//static void apply_parametric_invflow(float *y, float *x, int w, int h, int pd,
//		char *model_id, double *param, int nparam)
}

// paraflow A.PNG B.PNG model error
// model: traslation, euclidean, similar, affine, projective, {c,}parabolic, real1, real2
// error: l1, l2, linf, correlation, entropy
// cparabolic: radial distortion (x,y)->(x,y)+a*(x^2+y^2)*(x,y)
// parabolic: radial distortion with arbitrary center
// real1 = projective + parabolic
// real2 = parabolic + projective
int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s a b model \"pars\" error\n", *v);
		//                          0 1 2 3       4      5
		return EXIT_FAILURE;
	}
	int wx, hx, pdx, wy, hy, pdy;
	float *x = iio_read_image_float_vec(v[1], &wx, &hx, &pdx);
	float *y = iio_read_image_float_vec(v[2], &wy, &hy, &pdy);
	if (wx != wy || hx != hy || pdx != pdy)
		error("input size mismatch");

	int maxparam = 10;
	double param[maxparam];
	int nparams = parse_doubles(param, maxparam, v[4]);

	fprintf(stderr, "MODEL: %s\n", v[3]);
	FORI(nparams)
		fprintf(stderr, "PARAMETER[%d]: %g\n", i, param[i]);

	do_stuff(x, y, wx, hx, pdx, v[3], v[5], param, nparams);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
