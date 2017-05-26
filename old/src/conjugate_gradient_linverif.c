#include <stdlib.h>
#include "random.c"

static void vector_sum(double *xpy, double *x, double *y, int n)
{
	for (int i = 0; i < n; i++)
		xpy[i] = x[i] + y[i];
}

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
