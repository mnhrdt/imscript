#include <math.h>
#include "xmalloc.c"


typedef void (*linear_map_t)(double *y, double *x, int n, void *e);

static double scalar_product_rec(double *x, double *y, int n)
{
	return n > 0 ? *x * *y + scalar_product_rec(x+1, y+1, n-1) : 0;
}

static double scalar_product(double *x, double *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

#ifdef LINEAR_MAP_VERIFICATION
#include <stdlib.h>
#include "random.c"
static void fill_random_vector(double *x, int n)
{
	double mu = 1000*(random_uniform()-0.5);
	double sigma = 500*random_uniform()+0.1;
	for (int i = 0; i < n; i++)
		x[i] = mu + sigma*random_normal();
}

static void assert_symmetry(linear_map_t A, int n, void *e, int ntrials)
{
	int sn = n * sizeof(double);
	double *x = xmalloc(sn);
	double *y = xmalloc(sn);
	double *Ax = xmalloc(sn);
	double *Ay = xmalloc(sn);
	for (int i = 0; i < ntrials; i++) {
		fill_random_vector(x, n);
		fill_random_vector(y, n);
		A(Ax, x, n, e);
		A(Ay, y, n, e);
		double Axy = scalar_product(Ax, y, n);
		double xAy = scalar_product(x, Ay, n);
		double error = fabs(Axy-xAy);
		fprintf(stderr, "error symmetry = %g\n", error);
		if (error > 1e-12)
			fprintf(stderr, "WARNING: non symmetric!\n");
	}
	free(x);
	free(y);
	free(Ax);
	free(Ay);
}

static void assert_positivity(linear_map_t A, int n, void *e, int ntrials)
{
	int sn = n * sizeof(double);
	double *x = xmalloc(sn);
	double *Ax = xmalloc(sn);
	for (int i = 0; i < ntrials; i++) {
		fill_random_vector(x, n);
		A(Ax, x, n, e);
		double Axx = scalar_product(Ax, x, n);
		double pos = fabs(Axy-xAy);
		fprintf(stderr, "positivity = %g\n", pos);
		if (pos < -1e-12)
			fprintf(stderr, "WARNING: non positive!\n");
	}
	free(x);
	free(Ax);
}
#endif//LINEAR_MAP_VERIFICATION


#define FOR(i,n) for(int i = 0; i < n; i++)

void fancy_conjugate_gradient(double *x,
		linear_map_t A, double *b, int n, void *e,
		double *x0, int max_iter, double min_residual)
{
	double *r  = xmalloc(n * sizeof(double));
	double *p  = xmalloc(n * sizeof(double));
	double *Ap = xmalloc(n * sizeof(double));

	A(Ap, x, n, e);

	FOR(i,n) x[i] = x0[i];
	FOR(i,n) r[i] = b[i] - Ap[i];
	FOR(i,n) p[i] = r[i];

	for (int iter = 0; iter < max_iter; iter++) {
		A(Ap, p, n, e);
		double   App    = scalar_product(Ap, p, n);
		double   rr_old = scalar_product(r, r, n);
		double   alpha  = rr_old / App;
		FOR(i,n) x[i]   = x[i] + alpha * p[i];
		FOR(i,n) r[i]   = r[i] -alpha * Ap[i];
		double   rr_new = scalar_product(r, r, n);
		if (sqrt(rr_new) < min_residual)
			break;
		double   beta   = rr_new / rr_old;
		FOR(i,n) p[i]   = r[i] + beta * p[i];
	}

	free(r);
	free(p);
	free(Ap);
}


void conjugate_gradient(double *x, linear_map_t A, double *b, int n, void *e)
{
#ifdef LINEAR_MAP_VERIFICATION
	int number_of_trials = 10;
	assert_symmetry(A, n, e, number_of_trials);
	assert_positivity(A, n, e, number_of_trials);
#endif//LINEAR_MAP_VERIFICATION

	for (int i = 0; i < n; i++)
		x[i] = 0;
	int max_iter = n;
	double min_residual = 1e-12;

	fancy_conjugate_gradient(x, A, b, n, e, x, max_iter, min_residual);
}
