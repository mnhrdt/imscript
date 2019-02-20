// toolkit for dealing with RPC functions

#include <math.h> // sqrt

// just for debugging shit
#include <stdio.h>
static void matrix_print(double *A, int n, int m, char *name)
{
	double (*a)[m] = (void*)A;

	if (name)
		printf("%s =%s", name, n > 1 ? "\n" : "");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			printf("\t%g", a[i][j]);
		printf("\n");
	}
}


// cholesky decomposition (slow, trivial, only for small and full matrices)
// Note: it produces a lower triangular matrix B such that B*B'=A
static void cholesky(double *B, double *A, int n)
{
	double (*a)[n] = (void *)A;  // a[n][n]
	double (*b)[n] = (void *)B;  // a[n][n]

	for (int i = 0; i < n*n; i++)
		B[i] = 0;

	// standard Cholesky-Crout algorithm
	for (int i = 0; i < n; i++)
	for (int j = 0; j <= i; j++)
	{
		double r = 0;
		for (int k = 0; k < j; k++)
			r += b[i][k] * b[j][k];
		if (i == j)
			b[i][j] = sqrt(a[i][j] - r);
		else
			b[i][j] = ( a[i][j] - r ) / b[j][j];
	}

	// note: if the input is not a positive definite matrix,
	// this algorithm fills the last part of the matrix with NAN
}

// solve an upper triangular system A' * x = b
// where the matrix A is lower triangular
static void trisolve_upper(double *x, double *A, double *b, int n)
{
	double (*a)[n] = (void*)A;

	for (int k = n-1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int j = k+1; j < n; j++)
			x[k] -= a[j][k] * x[j];
		x[k] /= a[k][k];
	}
}

// solve a lower triangular system A * x = b
// where the matrix A is lower triangular
static void trisolve_lower(double *x, double *A, double *b, int n)
{
	double (*a)[n] = (void*)A;

	for (int k = 0; k < n; k++)
	{
		x[k] = b[k];
		for (int j = 0; j < k; j++)
			x[k] -= a[k][j] * x[j];
		x[k] /= a[k][k];
	}
}

// solve a symmetric, positive definite system
static void psymsolve(double *x, double *A, double *b, int n)
{
	// factor the matrix A = L * L'
	double L[n*n];
	cholesky(L, A, n);

	matrix_print(L, n, n, "L");

	// solve the two triangular systems L * L' * x = b
	double t[n];
	trisolve_lower(t, L, b, n);
	trisolve_upper(x, L, t, n);
}

static double scalar_product(double *x, double *y, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

// evaluate all 20 monomials of degree 3 in 3 variables
// (TODO: use a better conditioned basis)
static void eval_mon20(double m[20], double X[3])
{
	double x = X[0];
	double y = X[1];
	double z = X[2];
	double t[20] = {1, x, y, z, x*y,
		x*z, y*z, x*x, y*y, z*z,
		x*y*z, x*x*x, x*y*y, x*z*z, x*x*y,
		y*y*y, y*z*z, x*x*z, y*y*z, z*z*z};
	for (int i = 0; i < 20; i++)
		m[i] = t[i];
}

static void eval_mon4(double m[4], double x)
{
	m[0] = 1;
	m[1] = x;
	m[2] = x*x;
	m[3] = x*x*x;
}

// evaluate a polynomial of degree 3 in 3 variables
static double eval_pol20(double c[20], double X[3])
{
	double m[20];
	eval_mon20(m, X);
	return scalar_product(m, c, 20);
}

// evaluate a polynomial of degree 3 in 1 variable
static double eval_pol4(double c[4], double X)
{
	double m[4];
	eval_mon4(m, X);
	return scalar_product(m, c, 4);
}

double rpceval33(double p[20], double q[20], double x[3])
{
	return eval_pol20(p, x) / eval_pol20(q, x);
}

double rpceval13(double p[4], double q[4], double x)
{
	return eval_pol4(p, x) / eval_pol4(q, x);
}

double rpc33_error(double p[20], double q[20], double (*x)[3], double *f, int n)
{
	double v[n];
	for (int i = 0; i < n; i++)
		v[i] = rpceval33(p, q, x[i]);

	long double r = 0;
	for (int i = 0; i < n; i++)
		r = pow(f[i] - v[i], 2);
	return 100 * sqrt(r/n);
}

double rpc13_error(double p[4], double q[4], double *x, double *f, int n)
{
	double v[n];
	for (int i = 0; i < n; i++)
		v[i] = rpceval13(p, q, x[i]);

	long double r = 0;
	for (int i = 0; i < n; i++)
		r = pow(f[i] - v[i], 2);
	return 100 * sqrt(r/n);
}

// least-squares fit a 3rd degree RPC model to the given data values in 3D
void rpcfit33(double p[20], double q[20], double (*x)[3], double *f, int n)
{
	// Input: n data points/values (x,f)
	// Output: a rational function p/q of degree 3
	// Algorithm: minimize E(p,q) = w(x)|p(x)-q(x)f|^2
	//            iteratively reweighted by w=1/q^2
	//
	// Note: the minimization problem is linear, symmetric
	// and positive definite.  The solution is not 0 since
	// the zero-degree coefficient of q is fixed at 1.
	// Thus the system has size 39x39.
	//
	// Note: this algorithm fails in two cases
	// 1) when the data points fit exactly to a rational function of lower
	// degree (e.g., a pinhole camera)
	// 2) when there are not enough data points
	// In these cases the matrix of the linear problem is degenerate and
	// the current solver fails, giving all-NAN, parameters.
	//

	// initialize q
	p[0] = q[0] = 1;
	for (int i = 1; i < 20; i++)
		p[i] = q[i] = 0;

	printf("init err = %g\n", rpc33_error(p, q, x, f, n));

	// build matrix "M" of evaluated monomials
	double M[n][39];
	for (int i = 0; i < n; i++)
	{
		eval_mon20(M[i], x[i]);
		for (int j = 1; j < 20; j++)
			M[i][19+j] = -f[i] * M[i][j];
	}
	matrix_print((void*)M, n, 39, "M");

	int number_of_iterations = 13;
	for (int iteration = 0; iteration < number_of_iterations; iteration++)
	{
		// vector "w" of weights
		double w[n];
		for (int i = 0; i < n; i++)
			w[i] = pow( eval_pol20(q, x[i]), -2);
		matrix_print(w, 1, n, "W");

		// matrix "A" of the system
		double A[39][39] = {{0}};
		for (int i = 0; i < 39; i++)
		for (int j = 0; j < 39; j++)
		for (int k = 0; k < n ; k++)
			A[i][j] += w[k] * M[k][i] * M[k][j];
		matrix_print((void*)A, 39, 39, "A");

		// vector "b"
		double b[39] = {0};
		for (int i = 0; i < 39; i++)
		for (int j = 0; j < n ; j++)
			b[i] += w[j] * f[j] * M[j][i];
		matrix_print((void*)b, 1, 39, "b");

		// solve the linear problem
		double pq[39];
		psymsolve(pq, (void*)A, b, 39);
		matrix_print((void*)pq, 1, 39, "pq");

		// copy to output
		for (int i = 0; i < 20; i++)
			p[i] = pq[i];
		for (int i = 1; i < 20; i++)
			q[i] = pq[19+i];

		// debugging stuff
		printf("iter_%d err = %g\n", 1+iteration, rpc33_error(p, q, x, f, n));
		matrix_print(p, 1, 20, "p");
		matrix_print(q, 1, 20, "q");

	}
}


// least-squares fit a 3rd degree RPC model to the given data values in 1D
void rpcfit13(double p[4], double q[4], double *x, double *f, int n)
{
	// Input: n data points/values (x,f)
	// Output: a rational function p/q of degree 3
	// Algorithm: minimize E(p,q) = w(x)|p(x)-q(x)f|^2
	//            iteratively reweighted by w=1/q^2
	//
	// Note: the minimization problem is linear, symmetric
	// and positive definite.  The solution is not 0 since
	// the zero-degree coefficient of q is fixed at 1.
	// Thus the system has size 7x7.
	//

	// initialize q
	p[0] = q[0] = 1;
	for (int i = 1; i < 4; i++)
		p[i] = q[i] = 0;

	printf("init err = %g\n", rpc13_error(p, q, x, f, n));

	// build matrix "M" of evaluated monomials
	double M[n][7];
	for (int i = 0; i < n; i++)
	{
		eval_mon4(M[i], x[i]);
		for (int j = 1; j < 4; j++)
			M[i][3+j] = -f[i] * M[i][j];
	}
	matrix_print((void*)M, n, 7, "M");

	int number_of_iterations = 15;
	for (int iteration = 0; iteration < number_of_iterations; iteration++)
	{
		// vector "w" of weights
		double w[n];
		for (int i = 0; i < n; i++)
			w[i] = pow( eval_pol4(q, x[i]), -2);
		matrix_print(w, 1, n, "W");

		// matrix "A" of the system
		double A[7][7] = {{0}};
		for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
		for (int k = 0; k < n ; k++)
			A[i][j] += w[k] * M[k][i] * M[k][j];
		matrix_print((void*)A, 7, 7, "A");

		// vector "b"
		double b[7] = {0};
		for (int i = 0; i < 7; i++)
		for (int j = 0; j < n ; j++)
			b[i] += w[j] * f[j] * M[j][i];
		matrix_print((void*)b, 1, 7, "b");

		// solve the linear problem
		double pq[7];
		psymsolve(pq, (void*)A, b, 7);
		matrix_print((void*)pq, 1, 7, "pq");

		// copy to output
		for (int i = 0; i < 4; i++)
			p[i] = pq[i];
		for (int i = 1; i < 4; i++)
			q[i] = pq[3+i];

		printf("iter_%d err = %g\n", 1+iteration, rpc13_error(p, q, x, f, n));
		matrix_print(p, 1, 4, "p");
		matrix_print(q, 1, 4, "q");

	}
}


// tests

#ifdef MAIN_TEST

//  AB  =   A  *  B
// (np) = (nm) x (mp)
void matrix_product(double *AB, double *A, double *B, int n, int m, int p)
{
	double  (*a)[m] = (void*)A;  //  a[n][m]
	double  (*b)[p] = (void*)B;  //  b[m][p]
	double (*ab)[p] = (void*)AB; // ab[n][p]

	for (int i = 0; i < n*p; i++)
		AB[i] = 0;

	for (int i = 0; i < n; i++)
	for (int j = 0; j < p; j++)
	for (int k = 0; k < m; k++)
		ab[i][j] += a[i][k] * b[k][j];
}

void matrix_transpose(double *At, double *A, int n, int m)
{
	double (*at)[n] = (void*)At; // at[m][n]
	double  (*a)[m] = (void*)A;  //  a[n][m]

	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
		at[j][i] = a[i][j];
}

#include <assert.h>
#include "random.c"

static void test_solver(void)
{
	// test symmetric solver (solution is x = 1 2 3 )
	double A[9] = { 8, 2, 3, 2, 5, 0, 3, 0, 9 };
	double b[3] = {21, 12, 30};
	matrix_print(A, 3, 3, "A");
	matrix_print(b, 3, 1, "b");
	double x[3];
	psymsolve(x, A, b, 3);
	matrix_print(x, 3, 1, "x");
}

static void test_3d(void)
{
	// test RPC eval
	int n = 7; // number of samples along each dimension
	double X[n*n*n][3];
	int cx = 0;
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	for (int k = 0; k < n; k++)
	{
		X[cx][0] = -1 + i*2.0/(n-1);
		X[cx][1] = -1 + j*2.0/(n-1);
		X[cx][2] = -1 + k*2.0/(n-1);
		cx += 1;
	}
	assert(cx == n*n*n);

	//double P[20] = {0, 0.1, 0.2, 1, 0};
	//double Q[20] = {1, 0};
	double P[20];
	double Q[20];
	for (int i = 0; i < 20; i++)
	{
		P[i] = 0.0001 * random_normal();
		Q[i] = 0.0001 * random_normal();
	}
	P[0] = 0; P[1] = 0.1; P[2] = 0.2; P[3] = 1;
	Q[0] = 1;
	matrix_print(P, 1, 20, "P");
	matrix_print(Q, 1, 20, "Q");
	double f[n*n*n];
	for (int i = 0; i < n*n*n; i++)
		f[i] = rpceval33(P, Q, X[i]);
	for (int i = 0; i < n*n*n; i++)
		printf("X_%d =\t%g\t%g\t%g\t» %g\n",i,X[i][0],X[i][1],X[i][2],f[i]);

	double p[20] = {0};
	double q[20] = {0};
	rpcfit33(p, q, X, f, n*n*n);

	for (int i = 0; i < 20; i++) if (p[i]) printf("p[%d] = %g\n", i, p[i]);
	for (int i = 0; i < 20; i++) if (q[i]) printf("q[%d] = %g\n", i, q[i]);
}

static void test_1d(void)
{
	// test RPC eval
	int n = 67; // number of samples along each dimension
	double X[n];
	int cx = 0;
	for (int i = 0; i < n; i++)
		X[i] = -1 + i*2.0/(n-1);

	double P[4] = {0.01   , 1   , 0.012   , 0.023   };
	double Q[4] = {1   , 0.001   , 0.0002   , 0.0003   };
	//double P[4] = {0.00   , 1   , 0.000   , 0.000   };
	//double Q[4] = {1   , 0.000   , 0.0000   , 0.0000   };
	double f[n];
	for (int i = 0; i < n; i++)
		f[i] = rpceval13(P, Q, X[i]);
	for (int i = 0; i < n; i++)
		printf("X_%d =\t%g\t» %g\n",i,X[i],f[i]);

	double p[4] = {0};
	double q[4] = {0};
	rpcfit13(p, q, X, f, n);

	for (int i = 0; i < 4; i++) if (p[i]) printf("p[%d] = %g\n", i, p[i]);
	for (int i = 0; i < 4; i++) if (q[i]) printf("q[%d] = %g\n", i, q[i]);
}

int main(void)
{
	test_solver();
	test_3d();


	//double b[3] = {14, 23, 21};
	////double v[3] = {1, 10, 34};
	//double U[9] = { 1, 2, 3, 0, 4, 5, 0, 0, 7};
	//double L[9] = { 1, 0, 0, 0, 2, 0, 0, 0, 3};
	//matrix_transpose(L, U, 3, 3);
	//matrix_print(U, 3, 3, "U");
	//matrix_print(L, 3, 3, "L");

	//double x[3];
	//trisolve_upper(x, U, b, 3);
	//matrix_print(x, 3, 1, "x");
	//x[0] = x[1] = x[2] = NAN;
	//trisolve_lower(x, L, v, 3);
	//matrix_print(x, 3, 1, "x");
	//matrix_print(b, 3, 1, "b");
	//matrix_print(v, 3, 1, "v");

	//double C[9];
	//cholesky(C, A, 3);
	//matrix_print(C, 3, 3, "C");

	//double Ct[9];
	//matrix_transpose(Ct, C, 3, 3);
	//matrix_print(Ct, 3, 3, "Ct");

	//double CCt[9], CtC[9];
	//matrix_product(CCt, C, Ct, 3, 3, 3);
	//matrix_product(CtC, Ct, C, 3, 3, 3);
	//matrix_print(CCt, 3, 3, "CCt");
	//matrix_print(CtC, 3, 3, "CtC");

	//double t[3], tt[3];
	//trisolve_lower(t, C, b, 3);
	//matrix_print(t, 3, 1, "t");
	//trisolve_upper(tt, C, t, 3);
	//matrix_print(tt, 3, 1, "t");


	//double A[4] = { 1, 2, 3, 4 };
	////double B[8] = { 1,0, 1,0, 1,-1, 1,0 };
	////double A[4] = { 1, 0, 0, 0 };
	//double B[8] = { 0, 0, 0, 0,  1, 0, 0, 0 };
	//matrix_print(A, 1, 4, "A");
	//matrix_print(B, 4, 2, "B");
	//double AB[2];
	//matrix_product(AB, A, B, 1, 4, 2);
	//matrix_print(AB, 1, 2, "AB");

	////int n = 5;
//	double AB[9];
//	double BA[9];
//	matrix_transpose(B, A, 3, 3);
//	matrix_print(A, 3, 3, "A");
//	matrix_print(B, 3, 3, "B");
//	matrix_product(AB, A, B, 3, 3, 3);
//	matrix_product(BA, B, A, 3, 3, 3);
//	matrix_print(AB, 3, 3, "AB");
//	matrix_print(BA, 3, 3, "BA");
//	cholesky_old(C, AB, 3);
//	matrix_print(C, 3, 3, "C");
//	cholesky_old(C, BA, 3);
//	matrix_print(C, 3, 3, "C");
//
//	cholesky(C, AB, 3);
//	matrix_print(C, 3, 3, "C'");
//	cholesky(C, BA, 3);
//	matrix_print(C, 3, 3, "C'");

	return 0;
}



static void cholesky(double *L, double *A, int n);
static void lu(double *L, double *U, double *A, int n);
static void qr(double *Q, double *R, double *A, int n);
static void svd(double *U, double *S, double *V, double *A, int m, int n);


#endif//MAIN_TEST
