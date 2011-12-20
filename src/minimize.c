// A program for minimizing other programs
//
//
//
// 1. CONTEXT
//
// Suppose that you have a program (or a script file) that is called
// in the following way:
//
// 	./program opt1 opt2 ... optM par1 par2 ... parN
//
// where par1 ... parN are floating-point numbers and opt1 ... optN are
// arbitrary strings.
//
// Suppose that running this program prints a single floating point
// number to stdout, which is computed deterministically from its
// arguments.  This number is called the "result" of the program.
//
// Suppose that you want to find the set of values of par1 ... parN
// where the result attains its minimum (for a fixed choice of opt1
// ... optM).
//
// Then, this program is for you.  Although global minimization is
// unsolvable in general, doing it automatically is more confortable
// than doing it by hand, and faster than by exhaustive search on a
// fixed grid.
//
//
// 2. SYNTAX
//
// 	./minimize program "optset" "parset0" "parset1" ... "parsetN"
//
//

#define _POSIX_C_SOURCE 2

#include <stdio.h>
#include <stdlib.h>

#include "fail.c"



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

struct program {
	char *name;
	char *optset;
	int n;
};

static double run_program(struct program *p, double *x)
{
	int bufsize = 0x1000;
	char buf[bufsize];

	int n = snprintf(buf, bufsize, "%s %s", p->name, p->optset);
	for (int i = 0; i < p->n; i++)
		n += snprintf(buf+n, bufsize-n, " %f", x[i]);
	//snprintf(buf+n, bufsize-n, "");
	fprintf(stderr, "RUN %s\n", buf);

	FILE *f = popen(buf, "r");
	double y;
	if (1 != fscanf(f, "%lf\n", &y))
		fail("command [%s] did not print a number\n", buf);
	pclose(f);

	return y;
}

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

int main(int c, char *v[])
{
	if (c < 4) {
		fprintf(stderr, "usage:\n\t"
		"%s program \"optset\" \"parset0\" ... \"parsetN\"\n", *v);
	//       0  1         2          3               c-1
		return EXIT_FAILURE;
	}
	int nmax = c  - 3 - 1;
	struct program p = {.name = v[1], .optset = v[2], .n = nmax};
	for (int i = 3; i < c; i++) {
		double t[nmax];
		int n = parse_doubles(t, nmax, v[i]);
		if (n != nmax)
			fail("you must supply %d vectors of length %d,\n"
			"however, your %dth vector has length %d\n",
						nmax+1, nmax, i-2, n);
		double x = run_program(&p, t);
		fprintf(stderr, "returned x = %g\n", x);
	}

	return EXIT_SUCCESS;
}
