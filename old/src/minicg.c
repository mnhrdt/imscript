// a matrix-free linear solver


#include <math.h>   // only for "sqrt"
#include <stdlib.h> // only for "malloc"


// type of a function that computes a linear map on n-vectors
typedef void (*linear_map_t)(double *y, double *x, int n, void *e);

// compute the scalar product of two n-vectors
static double scalar_product(double *x, double *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

// macro for einstein notation
#define FOR(i,n) for(int i = 0; i < n; i++)

// general conjugate gradient method
void fancy_conjugate_gradient(double *x,
		linear_map_t A, double *b, int n, void *e,
		double *x0, int max_iter, double min_residual)
{
	double *r  = malloc(n * sizeof(double));
	double *p  = malloc(n * sizeof(double));
	double *Ap = malloc(n * sizeof(double));

	A(Ap, x0, n, e);

	FOR(i,n) x[i] = x0[i];
	FOR(i,n) r[i] = b[i] - Ap[i];
	FOR(i,n) p[i] = r[i];

	for (int iter = 0; iter < max_iter; iter++) {
		A(Ap, p, n, e);
		double   App    = scalar_product(Ap, p, n);
		double   rr_old = scalar_product(r, r, n);
		double   alpha  = rr_old / App;
		FOR(i,n) x[i]   = x[i] + alpha * p[i];
		FOR(i,n) r[i]   = r[i] - alpha * Ap[i];
		double   rr_new = scalar_product(r, r, n);
		fprintf(stderr, "n=%d iter=%d, rr_new=%g\n", n, iter, rr_new);
		if (sqrt(rr_new) < min_residual)
			break;
		double   beta   = rr_new / rr_old;
		FOR(i,n) p[i]   = r[i] + beta * p[i];
	}

	free(r);
	free(p);
	free(Ap);
}

#undef FOR

// conjugate gradient method with default parameters
void conjugate_gradient(double *x, linear_map_t A, double *b, int n, void *e)
{
	for (int i = 0; i < n; i++)
		x[i] = 0;
	int max_iter = n;
	double min_residual = 1e-10;

	fancy_conjugate_gradient(x, A, b, n, e, x, max_iter, min_residual);
}
