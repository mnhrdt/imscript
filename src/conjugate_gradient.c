// a matrix-free linear solver


#include <math.h>
#include "xmalloc.c"


typedef void (*linear_map_t)(double *y, double *x, int n, void *e);

static double scalar_product(double *x, double *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}




#define FOR(i,n) for(int i = 0; i < n; i++)

void fancy_conjugate_gradient(double *x,
		linear_map_t A, double *b, int n, void *e,
		double *x0, int max_iter, double min_residual)
{
	double *r  = xmalloc(n * sizeof(double));
	double *p  = xmalloc(n * sizeof(double));
	double *Ap = xmalloc(n * sizeof(double));

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
		fprintf(stderr, "iter=%d, rr_new=%g\n", iter, rr_new);
		if (sqrt(rr_new) < min_residual)
			break;
		double   beta   = rr_new / rr_old;
		FOR(i,n) p[i]   = r[i] + beta * p[i];
	}

	free(r);
	free(p);
	free(Ap);
}

#ifdef LINEAR_MAP_VERIFICATION
#include "conjugate_gradient_linverif.c"
#endif//LINEAR_MAP_VERIFICATION

void conjugate_gradient(double *x, linear_map_t A, double *b, int n, void *e)
{
#ifdef LINEAR_MAP_VERIFICATION
	int number_of_trials = 10;
	assert_linearity(A, n, e, number_of_trials);
	assert_symmetry(A, n, e, number_of_trials);
	assert_positivity(A, n, e, number_of_trials);
	if (n < 34) {
		print_matrices_stderr(A, b, n, e);
	}
#endif//LINEAR_MAP_VERIFICATION

	for (int i = 0; i < n; i++)
		x[i] = 0;
	int max_iter = n;
	double min_residual = 1e-10;

	fancy_conjugate_gradient(x, A, b, n, e, x, max_iter, min_residual);
}


#ifdef TEST_MAIN_CONJUGATE_GRADIENT

static void apply_matrix(double *y, double *x, int n, void *AA)
{
	double (*A)[n] = AA;
	for (int i = 0; i < n; i++)
	{
		y[i] = 0;
		for (int j = 0; j < n; j++)
			y[i] += x[j]*A[j][i];
	}
}

int main(void)
{
	double A[4] = {4, 2, 2, 3};
	double b[2] = {3,3};
	double x[2];
	conjugate_gradient(x, apply_matrix, b, 2, A);
	fprintf(stderr, "x = [%g %g]\n", x[0], x[1]);
	return 0;
}
#endif//TEST_MAIN_CONJUGATE_GRADIENT
