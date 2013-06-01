// a matrix-free linear solver


#include <math.h>
#include "xmalloc.c"


typedef void (*linear_map_t)(double *y, double *x, int n, void *e);

//static double scalar_product_rec(double *x, double *y, int n)
//{
//	return n > 0 ? *x * *y + scalar_product_rec(x+1, y+1, n-1) : 0;
//}

static double scalar_product(double *x, double *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

static void vector_sum(double *xpy, double *x, double *y, int n)
{
	for (int i = 0; i < n; i++)
		xpy[i] = x[i] + y[i];
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

static void assert_linearity(linear_map_t A, int n, void *e, int ntrials)
{
	int sn = n * sizeof(double);
	double *x = xmalloc(sn);
	double *y = xmalloc(sn);
	double *Ax = xmalloc(sn);
	double *Ay = xmalloc(sn);
	double *xpy = xmalloc(sn);
	double *AxpAy = xmalloc(sn);
	double *Axpy = xmalloc(sn);
	for (int i = 0; i < ntrials; i++) {
		fill_random_vector(x, n);
		fill_random_vector(y, n);
		vector_sum(xpy, x, y, n);
		A(Ax, x, n, e);
		A(Ay, y, n, e);
		A(Axpy, xpy, n, e);
		vector_sum(AxpAy, Ax, Ay, n);
		double error = 0;
		for (int k = 0; k < n; k++)
			error = hypot(error, fabs(AxpAy[k]-Axpy[k]));
		fprintf(stderr, "error linearity = %g\n", error);
		if (error > 1e-12)
			fprintf(stderr, "WARNING: non linear!\n");

	}
	free(x);
	free(y    );
	free(Ax   );
	free(Ay   );
	free(xpy  );
	free(AxpAy);
	free(Axpy );
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
		double pos = Axx;
		fprintf(stderr, "positivity = %g\n", pos);
		if (pos < -1e-12)
			fprintf(stderr, "WARNING: non positive!\n");
	}
	free(x);
	free(Ax);
}

static void print_matrices_stderr(linear_map_t A, double *b, int n, void *e)
{
	fprintf(stderr, "linear system of size %d\n", n);
	fprintf(stderr, "A =\n");
	for (int i = 0; i < n; i++)
	{
		double x[n], Ax[n];
		for (int j = 0; j < n; j++)
			x[j] = j == i;
		A(Ax, x, n, e);
		for (int j = 0; j < n; j++)
			fprintf(stderr, " %g", Ax[j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "b =\n");
	for (int j = 0; j < n; j++)
		fprintf(stderr, " %g", b[j]);
	fprintf(stderr, "\n");
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

#ifdef MAIN_CONJUGATE_GRADIENT

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
#endif//MAIN_CONJUGATE_GRADIENT
