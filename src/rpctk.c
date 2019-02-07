// toolkit for dealing with RPC functions

#include <stdio.h> // fprintf
#include <math.h> // sqrt

// conventions for matrix representation used in this file
//
// Matrices are represented as a pointer to a contiguous array of doubles,
// together with some ints that describe the dimensions and (optionally, the
// stride).

// cholesky decomposition (slow, trivial, only for small full matrices)
static void cholesky(double *B, double *A, int n)
{
	for (int i = 0; i < n*n; i++)
		B[i] = 0;

	// standard Cholesky-Crout algorithm
	for (int i = 0; i < n; i++)
	for (int j = 0; j <= i; j++)
	{
		double r = 0;
		for (int k = 0; k < j; k++)
			r += B[i*n + k] * B[j*n + k];
		if (i == j)
			B[i*n + j] = sqrt(A[i*n + j] - r);
		else
			B[i*n + j] = ( A[i*n + j] - r ) / B[j*n + j];
	}
}

// solve an upper triangular system A * x = b
// where the matrix A is upper triangular
static void trisolve_upper(double *x, double *A, double *b, int n)
{
	double (*a)[n] = (void*)A;

	for (int k = n-1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int j = k+1; j < n; j++)
			x[k] -= a[k][j] * x[j];
		x[k] /= a[k][k];
	}
}

// solve a lower triangular system A' * x = b
// where the matrix A is upper triangular
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
	// factor the matrix A = L' * L
	double L[n*n];
	cholesky(L, A, n);

	// solve the two triangular systems L' * L * x = b
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

// evaluate all 20 monomials of degree 3
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

// evaluate a polynomial of degree 3
static double eval_pol20(double c[20], double X[3])
{
	double m[20];
	eval_mon20(m, X);
	return scalar_product(m, c, 20);
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

	// initialize q
	q[0] = 1;
	for (int i = 1; i < 20; i++)
		q[i] = 0;

	// build matrix "M" of evaluated monomials
	double M[n][39];
	for (int i = 0; i < n; i++)
	{
		eval_mon20(M[i], x[i]);
		for (int j = 1; j < 20; j++)
			M[i][19+j] = -f[i] * M[i][j];
	}

	int number_of_iterations = 3;
	for (int iteration = 0; iteration < number_of_iterations; iteration++)
	{
		// vector "w" of weights
		double w[n];
		for (int i = 0; i < n; i++)
		{
			double r = eval_pol20(q, x[i]);
			w[i] = 1 / (r * r);
		}

		// matrix "A" of the system
		double A[39][39];
		for (int i = 0; i < 39; i++)
		for (int j = 0; j < 39; j++)
		{
			// TODO: compute only a triangle (A is symmetric)
			A[i][j] = 0;
			for (int k = 0; k < n; k++)
				A[i][j] += w[k] * M[k][i] * M[k][j];
		}

		// vector "b"
		double b[39];
		for (int i = 0; i < 39; i++)
		{
			b[i] = 0;
			for (int j = 0; j < n; j++)
				b[i] += w[j] * f[j] * M[j][i];
		}

		// solve the linear problem
		double pq[39];
		psymsolve(pq, (void*)A, b, 39);

		// assess the root mean square error at this iteration
		long double e = 0;
		for (int i = 0; i < n; i++)
		{
			double qq[20] = {1};
			for (int i = 0; i < 19; i++)
				qq[i + 1] = pq[i + 20];
			double px = eval_pol20(pq, x[i]);
			double qx = eval_pol20(qq, x[i]);
			double r = (px - f[i]*qx) / qx;
			e += r * r;
		}
		e = sqrt(e/n);
		fprintf(stderr, "err_%d = %g\n", iteration, (double)e);
	}
}


// tests

#ifdef MAIN_TEST

//  AB  =   A  *  B
// (np) = (nm) x (mp)
void matrix_product(double *AB, double *A, double *B, int n, int m, int p)
{
	double  (*a)[n] = (void*)A;  //  a[n][m]
	double  (*b)[m] = (void*)B;  //  b[m][p]
	double (*ab)[n] = (void*)AB; // ab[n][p]

	for (int i = 0; i < n*p; i++)
		AB[i] = 0;

	for (int i = 0; i < n; i++)
	for (int j = 0; j < p; j++)
	for (int k = 0; k < m; k++)
		ab[i][j] += a[i][k] * b[k][j];
}

void matrix_transpose(double *At, double *A, int n, int m)
{
	double (*at)[m] = (void*)At; // at[m][n]
	double  (*a)[n] = (void*)A;  //  a[n][m]

	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
		at[j][i] = a[i][j];
}

#include <stdio.h>
void matrix_print(double *A, int n, int m, char *name)
{
	double (*a)[n] = (void*)A;

	if (name)
		printf("%s =\n", name);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			printf("\t%g", a[j][i]);
		printf("\n");
	}
}

int main(void)
{
	double A[4] = { 1, 2, 3, 4 };
	//double B[8] = { 1,0, 1,0, 1,-1, 1,0 };
	//double A[4] = { 1, 0, 0, 0 };
	double B[8] = { 0, 0, 0, 0,  1, 0, 0, 0 };
	matrix_print(A, 1, 4, "A");
	matrix_print(B, 4, 2, "B");
	double AB[2];
	matrix_product(AB, A, B, 1, 4, 2);
	matrix_print(AB, 1, 2, "AB");

	////int n = 5;
	//double A[9] = {4, -1.01, -2,  0, 3, 10,  0, 0, 4};
	//double B[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	//double C[9];
	//double AB[9];
	//double BA[9];
	//matrix_transpose(B, A, 3, 3);
	//matrix_print(A, 3, 3, "A");
	//matrix_print(B, 3, 3, "B");
	//matrix_product(AB, A, B, 3, 3, 3);
	//matrix_product(BA, B, A, 3, 3, 3);
	//matrix_print(AB, 3, 3, "AB");
	//matrix_print(BA, 3, 3, "BA");
	//cholesky(C, AB, 3);
	//matrix_print(C, 3, 3, "C");
	//return 0;
}


#endif//MAIN_TEST
