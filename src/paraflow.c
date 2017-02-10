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

#include "fail.c"
#include "statistics.c"


#include "synflow_core.c"





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
	int n;
	char *model_id;
	char *error_id;
};

// this function is independent of GSL
static double eval_objective_function(struct problem_data *p, double *v, int nv)
{
	float *x = p->x;
	float *y = p->y;
	int w = p->w;
	int h = p->h;
	int pd = p->pd;
	float *bp = xmalloc(w * h * pd * sizeof*bp);
	struct flow_model fm[1];
	produce_flow_model(fm, v, nv, p->model_id, w, h);
	transform_back(bp, fm, y, w, h, pd);
	//apply_parametric_invflow(bp, y, w, h, pd, p->model_id, v, nv);
	//iio_save_image_float_vec("/tmp/merdota.tiff", bp, w, h, pd);
	double r = evaluate_error_between_images(bp, x, w, h, pd, p->error_id);
	free(bp);
	return r;
}

// this function is called by the GSL minimzator
static double objective_function(const gsl_vector *v, void *pp)
{
	struct problem_data *p = pp;
	assert(p->n == (int)v->size);
	double point[v->size];
	FORI((int)v->size)
		point[i] = gsl_vector_get(v, i);
	double r = eval_objective_function(p, point, v->size);
	return r;
}

//SMART_PARAMETER(SIMPLEX_NEQ,5)
#define SIMPLEX_NEQ 5

static bool are_equal(double a, double b)
{
	return fabs(b-a) < 0.000000000001;
}

// external interface to GSL code (for 4D vectors)
static int minimize_objective_function(void *pp,
		float *result, float *starting_point, float *stepsize)
{
	struct problem_data *p = pp;
	int dim = p->n;

	gsl_vector *ss = gsl_vector_alloc(dim);
	gsl_vector *x = gsl_vector_alloc(dim);

	FORI(dim) gsl_vector_set(ss, i, stepsize[i]);
	FORI(dim) gsl_vector_set(x, i, starting_point[i]);

	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, dim);

	gsl_multimin_function f = {.f=objective_function, .n=dim, .params=pp};

	gsl_multimin_fminimizer_set(s, &f, x, ss);

	int status, iter = 0;
	int scount = 0; // counts the number of repeated values
	double oldval = INFINITY;
	do {
		iter += 1;

		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;

		double size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 0.00001);

		if (status == GSL_SUCCESS)
			fprintf(stderr, "converged to minimum!\n");

		fprintf(stderr, "iter = %d (", iter);
		FORI(dim) fprintf(stderr, "%g ", gsl_vector_get(s->x, i));
		fprintf(stderr, ") f=%g size=%g\n", s->fval, size);

		//scount = are_equal(oldval, s->fval) ? scount+1 : 0;
		//if (scount > SIMPLEX_NEQ) break;
		oldval = s->fval;

	} while (status == GSL_CONTINUE && iter < 490);

	FORI(dim) result[i] = gsl_vector_get(s->x, i);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	return status;
}

//
// END OF GSL STUFF
//
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

static double corrf(float *x, float *y, int n)
{
	double mx = 0; FORI(n) mx += x[i]; mx /= n;
	double my = 0; FORI(n) my += y[i]; my /= n;
	double sx = 0; FORI(n) sx = hypot(sx, x[i] - mx); sx /= sqrt(n-1);
	double sy = 0; FORI(n) sy = hypot(sy, y[i] - my); sy /= sqrt(n-1);
	double cc = 0; FORI(n) cc += (x[i] - mx)*(y[i] - my);
	cc /= (n-1)*sx*sy;
	return cc;
}

static double kendall(float *x, float *y, int n)
{
	long double concordant = 0;
	long double discordant = 0;
	long double concordantp = 0;
	long double discordantp = 0;
	long long int np = 0;
	for (int i = 0; i < n; i++)
	for (int j = 0; j < i; j++)
	{
		if (x[i] > x[j] && y[i] > y[j]) concordant += 1;
		if (x[i] < x[j] && y[i] < y[j]) concordantp += 1;
		if (x[i] > x[j] && y[i] < y[j]) discordant += 1;
		if (x[i] < x[j] && y[i] > y[j]) discordantp += 1;
		np += 1;
	}
	long double npairs = n;
	npairs = 0.5 * npairs * (npairs-1);
	long double nonties = concordant+concordantp+discordant+discordantp;
	fprintf(stderr, "n = %d\n", n);
	fprintf(stderr, "np = %lld\n", np);
	fprintf(stderr, "npairs = %Lf\n", npairs);
	fprintf(stderr, "concordant = %Lf\n", concordant);
	fprintf(stderr, "discordant = %Lf\n", discordant);
	fprintf(stderr, "concordantp = %Lf\n", concordantp);
	fprintf(stderr, "discordantp = %Lf\n", discordantp);
	fprintf(stderr, "real concordant = %Lf\n", concordant+concordantp);
	fprintf(stderr, "real discordant = %Lf\n", discordant+discordantp);
	fprintf(stderr, "nonties = %Lf\n", nonties);
	//double kk = (concordant +concordantp - discordant - discordantp) / npairs;
	double kk = (concordant +concordantp - discordant - discordantp) / nonties;
	fprintf(stderr, "kk = %lf\n", kk);
	return kk;
}

#define ERRPCWINRADIUS 7

static double error_l2_at_pcpoint(float *xx, float *yy, int w, int h, int pd,
		float p[2])
{
	float (*x)[w][pd] = (void*)xx;
	float (*y)[w][pd] = (void*)yy;

	int i0 = w*p[0]/100;
	int j0 = h*p[1]/100;
	int i_start = i0 - ERRPCWINRADIUS;
	int i_end   = i0 + ERRPCWINRADIUS;
	int j_start = j0 - ERRPCWINRADIUS;
	int j_end   = j0 + ERRPCWINRADIUS;

	if (i_start < 0 || j_start < 0 || i_end >= w || j_end >= h)
		return -1;

	double r = 0;
	for (int j = j_start; j <= j_end; j++)
	for (int i = i_start; i <= i_end; i++)
	{
		double t = 0;
		FORL(pd)
			t = hypot(t, x[j][i][l] - y[j][i][l]);
		r += t*t;
	}

	return r;
}

static double error_l2_at_point(float *xx, float *yy, int w, int h, int pd,
		int p[2])
{
	float (*x)[w][pd] = (void*)xx;
	float (*y)[w][pd] = (void*)yy;


	int i0 = p[0];
	int j0 = p[1];
	int i_start = i0 - ERRPCWINRADIUS;
	int i_end   = i0 + ERRPCWINRADIUS;
	int j_start = j0 - ERRPCWINRADIUS;
	int j_end   = j0 + ERRPCWINRADIUS;

	if (i_start < 0 || j_start < 0 || i_end >= w || j_end >= h)
		return -1;

	double r = 0;
	for (int j = j_start; j <= j_end; j++)
	for (int i = i_start; i <= i_end; i++)
	{
		double t = 0;
		FORL(pd)
			t = hypot(t, x[j][i][l] - y[j][i][l]);
		r += t*t;
	}

	return r;
}

static double error_corr_at_point(float *xx, float *yy, int w, int h, int pd,
		int p[2])
{
	float (*x)[w][pd] = (void*)xx;
	float (*y)[w][pd] = (void*)yy;


	int i0 = p[0];
	int j0 = p[1];
	int i_start = i0 - ERRPCWINRADIUS;
	int i_end   = i0 + ERRPCWINRADIUS;
	int j_start = j0 - ERRPCWINRADIUS;
	int j_end   = j0 + ERRPCWINRADIUS;

	if (i_start < 0 || j_start < 0 || i_end >= w || j_end >= h)
		return -1;

	double r = 0;
	float t[pd*(2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1)];
	float s[pd*(2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1)];
	int n = 0;
	for (int j = j_start; j <= j_end; j++)
	for (int i = i_start; i <= i_end; i++)
	{
		FORL(pd) {
			t[n] = x[j][i][l];
			s[n] = y[j][i][l];
			n += 1;
		}
	}
	assert(n == pd * (2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1));

	double k = corrf(t, s, n);
	assert(k >= 0);
	assert(k <= 1);
	return 1 - k;
}

static double error_kendall_at_pcpoint(float *xx, float *yy,
		int w, int h, int pd,
		float p[2])
{
	float (*x)[w][pd] = (void*)xx;
	float (*y)[w][pd] = (void*)yy;

	int i0 = w*p[0]/100;
	int j0 = h*p[1]/100;
	int i_start = i0 - ERRPCWINRADIUS;
	int i_end   = i0 + ERRPCWINRADIUS;
	int j_start = j0 - ERRPCWINRADIUS;
	int j_end   = j0 + ERRPCWINRADIUS;

	if (i_start < 0 || j_start < 0 || i_end >= w || j_end >= h)
		return -1;

	double r = 0;
	float t[pd*(2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1)];
	float s[pd*(2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1)];
	int n = 0;
	for (int j = j_start; j <= j_end; j++)
	for (int i = i_start; i <= i_end; i++)
	{
		FORL(pd) {
			t[n] = x[j][i][l];
			s[n] = y[j][i][l];
			n += 1;
		}
	}
	assert(n == pd * (2*ERRPCWINRADIUS+1)*(2*ERRPCWINRADIUS+1));

	double k = kendall(t, s, n);
	assert(k >= 0);
	assert(k <= 1);
	return 1 - k;
}

static double error_at_pcpoints(float *xx, float *yy, int w, int h, int pd,
		float (*pcpoints)[2], int n,
		double (f)(float*,float*,int,int,int,float*) )
{
	double r = 0;
	int m = 0;
	FORI(n) {
		double t = f(xx, yy, w, h, pd, pcpoints[i]);
		if (t >= 0) {
			r += t;
			m += 1;
		}
	}
	if (!m) fail("could not evaluate error at any of the given points!");
	return r/m;
}

static double error_at_points(float *xx, float *yy, int w, int h, int pd,
		int (*points)[2], int n,
		double (f)(float*,float*,int,int,int,int*) )
{
	double r = 0;
	int m = 0;
	FORI(n) {
		double t = f(xx, yy, w, h, pd, points[i]);
		if (t >= 0) {
			r += t;
			m += 1;
		}
	}
	if (m != 8) fail("m != 8 (%d)", m);
	if (!m) fail("could not evaluate error at any of the given points!");
	return r/m;
}

static double evaluate_error_between_images(float *xx, float *yy,
		int w, int h, int pd, char *error_id)
{
	if (0 == strcmp(error_id, "l2points")) {
		int points[8][2] = {
			{169,978},
			{127,900},
			{469,943},
			{533,1046},
			{169,157},
			{393,241},
			{411,104},
			{161,395}
		};
		double r = error_at_points(xx, yy, w, h, pd, points, 8,
				error_l2_at_point);
		return r;
	}
	if (0 == strcmp(error_id, "corrpoints")) {
		int points[7][2] = {
			{169,978},
			{127,900},
			{469,943},
			{533,1046},
			{169,157},
			{393,241},
			{411,104}
		};
		double r = error_at_points(xx, yy, w, h, pd, points, 7,
				error_corr_at_point);
		return r;
	}
	if (0 == strcmp(error_id, "l2pcpoints")) {
		float pcpoints[4][2] = {{20,20},{20,80},{80,20},{80,80}};
		double r = error_at_pcpoints(xx, yy, w, h, pd, pcpoints, 4,
				error_l2_at_pcpoint);
		return r;
	}
	if (0 == strcmp(error_id, "kendallpcpoints")) {
		float pcpoints[4][2] = {{20,20},{20,80},{80,20},{80,80}};
		double r = error_at_pcpoints(xx, yy, w, h, pd, pcpoints, 4,
				error_kendall_at_pcpoint);
		return r;
	}
	//if (0 == strcmp(error_id, "corrpoints")) {
	//	float pcpoints[4][2] = {{20,20},{20,80},{80,20},{80,80}};
	//	double r = error_corr_at_pcpoints(x, y, w, h, pd, pcpoints, 4);
	//	return r;
	//}

	float (*x)[w][pd] = (void*)xx;
	float (*y)[w][pd] = (void*)yy;


	int passepartout = 8;
	double r = 0;
	int n = w*h*pd;
	float (*d)[w][pd] = xmalloc(n*sizeof(float));
	int nvals = 0;
	float *xvals = xmalloc(n*sizeof*xvals);
	float *yvals = xmalloc(n*sizeof*yvals);
	FORI(n) d[0][0][i] = 0;//fabs(x[i] - y[i]);
	for (int j = passepartout; j < h-passepartout; j++)
	for (int i = passepartout; i < w-passepartout; i++)
		FORL(pd) {
			d[j][i][l] = fabs(x[j][i][l] - y[j][i][l]);
			xvals[nvals] = x[j][i][l];
			yvals[nvals] = y[j][i][l];
			if (xvals[nvals] && yvals[nvals])
				nvals += 1;
		}

	//iio_save_image_float_vec("/tmp/diff.tiff", d[0][0], w, h, pd);
	if (0 == strcmp(error_id, "l2")) {
		r = 0;
		FORI(n)
			r = hypot(r, d[0][0][i]);
	} else if (0 == strcmp(error_id, "l1")) {
		r = 0;
		FORI(n)
			r = r + d[0][0][i];
	} else if (0 == strcmp(error_id, "linf")) {
		FORI(n)
			if (d[0][0][i] > r)
				r = d[0][0][i];
	} else if (0 == strcmp(error_id, "corr")) {
		r = corrf(xvals, yvals, nvals);
		r = fabs(r);
	} else if (0 == strcmp(error_id, "kendall")) {
		r = kendall(xvals, yvals, nvals);
		r = fabs(r);
		fail("kendall computation too slow!");
	} else fail("bad error id \"%s\"", error_id);
	//struct statistics_float s;
	//statistics_getf(&s, d, n);
	//print_stats(stderr, &s, "fabs(diff)");
	free(xvals);
	free(yvals);
	free(d);
	return r;
}


static void do_stuff(float *x, float *y, int w, int h, int pd,
		char *model_id, char *error_id,
		double *param, int nparams)
{
	struct problem_data p[1];
	p->w = w;
	p->h = h;
	p->pd = pd;
	p->x = x;
	p->y = y;
	p->n = nparams;
	p->model_id = model_id;
	p->error_id = error_id;

//	double r = eval_objective_function(p, param, nparams);
//static int minimize_objective_function(void *pp,
//		float *result, float *starting_point, float *stepsize)
	float result[nparams];
	float starting_point[nparams];
	float stepsize[nparams];
	FORI(nparams) starting_point[i] = param[i];
	FORI(nparams) stepsize[i] = 1;
	if (0 == strcmp(model_id, "affine")) {
		assert(nparams == 6);
		stepsize[0] = 0.0001;
		stepsize[1] = 0.0001;
		stepsize[2] = 0.1;
		stepsize[3] = 0.0001;
		stepsize[4] = 0.0001;
		stepsize[5] = 0.1;
	}
	int r= minimize_objective_function(p, result, starting_point, stepsize);

	fprintf(stderr, "gsl exit status = %d\n", r);
	FORI(nparams) fprintf(stderr, "result[%d] = %g\n", i, result[i]);

	{
		struct flow_model fm[1];
		double rrr[nparams]; FORI(nparams) rrr[i]=result[i];
		produce_flow_model(fm, rrr, nparams, model_id, w, h);
		float *flo = xmalloc(w*h*2*sizeof*flo);
		fill_flow_field(flo, fm, w, h);
		iio_save_image_float_vec("/tmp/fff.tiff", flo, w, h, 2);
		free(flo);
	}

	//fprintf(stderr, "r = %lf\n", r);
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
		fail("input size mismatch");

	int maxparam = 40;
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
