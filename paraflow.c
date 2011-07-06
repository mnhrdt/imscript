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
#include "statistics.c"



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


static void affine_map(double y[2], double A[6], double x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

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
			float a = x[iq][ip][l];
			float b = x[iq+1][ip][l];
			float c = x[iq][ip+1][l];
			float d = x[iq+1][ip+1][l];
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			result[l] = v;
		}
	}
}

static void apply_affinity(float *py, float *x, int w, int h, int pd,
		double A[6])
{
	float (*y)[w][pd] = (void*)py;
	FORJ(h) FORI(w) {
		double p[2] = {i, j}, q[2];
		affine_map(q, A, p);
		float val[pd];
		general_interpolate(val, x, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

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
	FORJ(h) FORI(w) {
		double p[2] = {i, j}, q[2];
		projective_map(q, H, p);
		float val[pd];
		general_interpolate(val, x, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}



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
	} else if (0 == strcmp(model_id, "traslation_s")) {
		assert(nparam == 2);
		double a[6] = {1, 0, w*param[0], 0, 1, h*param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "euclidean")) {
		assert(nparam == 3);
		double a[6] = {cos(param[2]), sin(param[2]), param[0],
			-sin(param[2]), cos(param[2]), param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "euclidean_s")) {
		assert(nparam == 3);
		double a[6] = {cos(param[2]), sin(param[2]), w*param[0],
			-sin(param[2]), cos(param[2]), h*param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "similar")) {
		assert(nparam == 4);
		double r = param[3];
		double a[6] = {r*cos(param[2]), r*sin(param[2]), param[0],
			-r*sin(param[2]), r*cos(param[2]), param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "similar_s")) {
		assert(nparam == 4);
		double r = param[3]/w;
		double a[6] = {r*cos(param[2]), r*sin(param[2]), w*param[0],
			-r*sin(param[2]), r*cos(param[2]), h*param[1]};
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "affine")) {
		assert(nparam == 6);
		double a[6]; FORI(6) a[i] = param[i];
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "affine_s")) {
		error("not implemented");
		assert(nparam == 6);
		double a[6];
		apply_affinity(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "projective")) {
		assert(nparam == 8);
		double a[9]; FORI(8) a[i] = param[i];
		a[8] = 1;
		apply_homography(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "projective_s")) {
		error("not implemented");
		assert(nparam == 8);
		double a[9];
		a[8] = 1;
		apply_homography(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "cparabolic")) {
		assert(nparam == 1);
		double a[3] = {w/2.0, h/2.0, param[0]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "cparabolic_s")) {
		assert(nparam == 1);
		double a[3] = {w/2.0, h/2.0, w*param[0]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "parabolic")) {
		assert(nparam == 3);
		double a[3] = {param[0], param[1], param[2]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "parabolic_s")) {
		assert(nparam == 3);
		double a[3] = {w*param[0], h*param[1], w*param[2]};
		apply_radialpol(y, x, w, h, pd, a);
	} else if (0 == strcmp(model_id, "real1")) {
		error("not yet implemented");
	} else if (0 == strcmp(model_id, "real2")) {
		error("not yet implemented");
	} else error("unrecognized flow model %s", model_id);
}

static double evaluate_error_between_images(float *x, float *y,
		int w, int h, int pd, char *error_id);

//
// HERE BE DRAGONS (GSL multidimensional minimization)
//

struct problem_data {
	// image pair
	int w, h, pd;
	float *x, *y;

	// model metadata
	char *model_id;
	char *error_id;
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
	//apply_parametric_invflow(iy, y, w, h, pd, model_id, param, nparams);
	iio_save_image_float_vec("/tmp/merdota.tiff", bp, p->w, p->h, p->pd);
	double r =
       	evaluate_error_between_images(bp, p->x, p->w, p->h, p->pd, p->error_id);
	free(bp);
	return r;
}

//SMART_PARAMETER(SIMPLEX_NEQ,5)

static bool are_equal(double a, double b)
{
	return fabs(b-a) < 0.000000000001;
}

// external interface to GSL code (for 4D vectors)
//static int minimize_objective_function(float *vmin, float *vzero, void *p)
//{
//	int dim = 4;
//	const gsl_multimin_fminimizer_type *T =
//		gsl_multimin_fminimizer_nmsimplex;
//	gsl_multimin_fminimizer *s = NULL;
//	gsl_multimin_function f = {.f=objective_function, .n=dim, .params=p};
//
//	gsl_vector *ss = gsl_vector_alloc(dim);
//	gsl_vector *x = gsl_vector_alloc(dim);
//
//	gsl_vector_set(ss, YOU_ANGLE, 1.0);
//	gsl_vector_set(ss, YOU_DIST, 0.05);
//	gsl_vector_set(ss, YOU_T, 0.05);
//	gsl_vector_set(ss, YOU_H, 0.05);
//
//	FORI(dim) gsl_vector_set(x, i, vzero[i]);
//
//	s = gsl_multimin_fminimizer_alloc(T, dim);
//	gsl_multimin_fminimizer_set(s, &f, x, ss);
//
//	int status, iter = 0;
//	int scount = 0; // counts the number of repeated values
//	double oldval = INFINITY;
//	do {
//		iter += 1;
//
//		status = gsl_multimin_fminimizer_iterate(s);
//		if (status) break;
//
//		double size = gsl_multimin_fminimizer_size(s);
//		status = gsl_multimin_test_size(size, 0.001);
//
//		if (status == GSL_SUCCESS)
//			fprintf(stderr, "converged to minimum!\n");
//
//		fprintf(stderr, "iter = %d (", iter);
//		FORI(dim) fprintf(stderr, "%g ", gsl_vector_get(s->x, i));
//		fprintf(stderr, ") f=%g size=%g\n", s->fval, size);
//
//		scount = are_equal(oldval, s->fval) ? scount+1 : 0;
//		if (scount > SIMPLEX_NEQ()) break;
//		oldval = s->fval;
//
//	} while (status == GSL_CONTINUE && iter < 50);
//
//	FORI(dim) vmin[i] = gsl_vector_get(s->x, i);
//
//	gsl_vector_free(x);
//	gsl_vector_free(ss);
//	gsl_multimin_fminimizer_free(s);
//
//	return status;
//}
//
////
//// END OF GSL STUFF
////
//
//
//
//SMART_PARAMETER_INT(MULTISCALE_START,5)
//
//#define NSCALES 8 //should depend on the size of the frames
//typedef struct {
//	imatge t[NSCALES];
//} multiscale;
//
//static void find_grass(imatge_rgba *outie, imatge *grass)
//{
//	multiscale m[1];
//	compute_median_multiscales(m, grass);
//
//	you u[1]={{.v={0,0,0,0},.corner=false}};
//	find_best_you_exhaustively(u, m->t + MULTISCALE_START());
//
//	float bestu[4]; FORI(4) bestu[i] = u->v[i];
//
//	minimize_objective_function(bestu, bestu, m->t + 4);
//	minimize_objective_function(bestu, bestu, m->t + 3);
//	minimize_objective_function(bestu, bestu, m->t + 2);
//	minimize_objective_function(bestu, bestu, m->t + 1);
//	minimize_objective_function(bestu, bestu, m->t + 0);
//
//	FORI(4) u->v[i] = bestu[i];
//
//	//you un[1];
//	//quadratic_amelioration(un, u, m->t + 6);
//	//{
//	//	float l[2][3];
//	//	get_you_cornerlines(l[0], l[1], grass, u, true);
//	//}
//
//	colorize_u(outie, grass, u);
//
//	FORI(NSCALES)
//		allibera_imatge(m->t + i);
//}

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

static double evaluate_error_between_images(float *xx, float *yy,
		int w, int h, int pd, char *error_id)
{
	float (*x)[w] = (void*)xx;
	float (*y)[w] = (void*)yy;
	// TODO: omit a passepartout
	int passepartout = 2;
	double r;
	int n = w*h*pd;
	float (*d)[w] = xmalloc(n*sizeof*d);
	FORI(n) d[0][i] = 0;//fabs(x[i] - y[i]);
	for (int j = passepartout; j < h-passepartout; j++)
	for (int i = passepartout; i < w-passepartout; i++)
		d[j][i] = fabs(x[j][i] - y[j][i]);

	iio_save_image_float_vec("/tmp/diff.tiff", d[0][0], w, h, pd);
	if (0 == strcmp(error_id, "l2")) {
		r = 0;
		FORI(n)
			r = hypot(r, d[0][i]);
	} else error("bad error id \"%s\"", error_id);
	//struct statistics_float s;
	//statistics_getf(&s, d, n);
	//print_stats(stderr, &s, "fabs(diff)");
	free(d);
	return r;
}


static void do_stuff(float *x, float *y, int w, int h, int pd,
		char *model_id, char *error_id,
		double *param, int nparams)
{
	float *iy = xmalloc(w*h*pd*sizeof*iy);
	apply_parametric_invflow(iy, y, w, h, pd, model_id, param, nparams);
	iio_save_image_float_vec("/tmp/merdota.tiff", iy, w, h, pd);
	double r = evaluate_error_between_images(iy, x, w, h, pd, error_id);
	free(iy);
	fprintf(stderr, "r = %lf\n", r);
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
