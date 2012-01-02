#ifndef _MINIMIZE_GSL
#define _MINIMIZE_GSL

// API
typedef double (objective_function)(double *x, int n, void *data);

// API
int minimize_objective_function(double *result, double *first, double *step,
		objective_function *f, int n, void *data);


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <gsl/gsl_multimin.h>



#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif



#include "smapa.h"

SMART_PARAMETER(GSL_MIN_MAXITER,500)
SMART_PARAMETER(GSL_MIN_DELTA,0.0001)
SMART_PARAMETER_SILENT(GSL_MIN_SIMPLEX_NEQ,5)
SMART_PARAMETER_SILENT(GSL_MIN_EPSILON,0.000000001)
//#define GSL_MIN_SIMPLEX_NEQ 5
//#define GSL_MIN_MAXITER 500
//#define GSL_MIN_DELTA 0.0001
//#define GSL_MIN_EPSILON 0.000000000001

static bool are_equal(double a, double b)
{
	return fabs(b-a) < GSL_MIN_EPSILON();
}


struct enveloped_function {
	objective_function *f;
	int n;
	void *data;
};

static double envelop_objective_function(const gsl_vector *v, void *pp)
{
	struct enveloped_function *p = pp;
	assert(p->n == (int)v->size);
	int n = p->n;
	double point[n];
	FORI(n)
		point[i] = gsl_vector_get(v, i);
	double r = (p->f)(point, n, p->data);
	return r;
}

// external interface to GSL code
int minimize_objective_function(double *result, double *first, double *step,
		objective_function *f, int n, void *data)
{
	gsl_vector *ss = gsl_vector_alloc(n);
	gsl_vector *x = gsl_vector_alloc(n);

	FORI(n) gsl_vector_set(ss, i, step[i]);
	FORI(n) gsl_vector_set(x, i, first[i]);

	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, n);

	struct enveloped_function e[1] = {{
		.f = f,
		.n = n,
		.data = data }};

	gsl_multimin_function ff = {
		.f = envelop_objective_function,
		.n = n,
		.params = e};

	gsl_multimin_fminimizer_set(s, &ff, x, ss);

	int status, iter = 0;
	int scount = 0; // counts the number of repeated values
	double oldval = INFINITY;
	do {
		iter += 1;

		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;

		double size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, GSL_MIN_DELTA());

		if (status == GSL_SUCCESS)
			fprintf(stderr, "converged to minimum!\n");

		fprintf(stderr, "iter = %d (", iter);
		FORI(n) fprintf(stderr, "%g ", gsl_vector_get(s->x, i));
		fprintf(stderr, ") f=%g size=%g\n", s->fval, size);

		if (GSL_MIN_SIMPLEX_NEQ()>0) {
			scount = are_equal(oldval, s->fval) ? scount+1 : 0;
			if (scount > GSL_MIN_SIMPLEX_NEQ()) break;
		}
		oldval = s->fval;

	} while (status == GSL_CONTINUE && iter < GSL_MIN_MAXITER());

	FORI(n) result[i] = gsl_vector_get(s->x, i);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	return status;
}
#endif//_MINIMIZE_GSL
