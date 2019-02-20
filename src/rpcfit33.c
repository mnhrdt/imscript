
#include <math.h> // sqrtl

// cholesky decomposition (slow, trivial, only for small and full matrices)
// Note: it produces a lower triangular matrix B such that B*B'=A
static void cholesky(long double *B, long double *A, int n)
{
	long double (*a)[n] = (void *)A;  // a[n][n]
	long double (*b)[n] = (void *)B;  // a[n][n]

	for (int i = 0; i < n*n; i++)
		B[i] = 0;

	// standard Cholesky-Crout algorithm
	for (int i = 0; i < n; i++)
	for (int j = 0; j <= i; j++)
	{
		long double r = 0;
		for (int k = 0; k < j; k++)
			r += b[i][k] * b[j][k];
		if (i == j)
			b[i][j] = sqrtl(a[i][j] - r);
		else
			b[i][j] = ( a[i][j] - r ) / b[j][j];
	}

	// NOTE: if the input is not a positive definite matrix,
	// this algorithm fills the last part of the matrix with NAN
	//
	// NOTE2: the Cholesky factorization is used below to compute a
	// Moore-Penrose pseudo-inverse.  In that context, using a full SVD is
	// a bit of overkill (although correct).  It seems more natural to use
	// a modified version of the Cholesky-Crout algorithm that detects when
	// "r" reaches zero and acts accordingly (providing a non-unique but
	// useful factorization).  See for example the algorithm here:
	// P. Corrieu, "Fast Computation of Moore-Penrose Pseudoinverses",
	// Neural Information Processing -- Letters and Reviews 2005
	// (which is a slight variation of the present function).
}

// solve an upper triangular system A' * x = b
// where the matrix A is lower triangular
static void solve_upper(long double *x, long double *A, long double *b, int n)
{
	long double (*a)[n] = (void*)A;

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
static void solve_lower(long double *x, long double *A, long double *b, int n)
{
	long double (*a)[n] = (void*)A;

	for (int k = 0; k < n; k++)
	{
		x[k] = b[k];
		for (int j = 0; j < k; j++)
			x[k] -= a[k][j] * x[j];
		x[k] /= a[k][k];
	}
}

// solve a symmetric, positive definite system
static void solve_spd(long double *x, long double *A, long double *b, int n)
{
	// factor the matrix A = L * L'
	long double L[n*n];
	cholesky(L, A, n);

	// solve the two triangular systems L * L' * x = b
	long double t[n];
	solve_lower(t, L, b, n);
	solve_upper(x, L, t, n);
}

static long double scalar_product(long double *x, long double *y, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

// evaluate all 20 monomials of degree 3 in 3 variables
// (TODO: use a better conditioned basis)
static void eval_mon20(long double m[20], long double X[3])
{
	long double x = X[0];
	long double y = X[1];
	long double z = X[2];
	long double t[20] = {1, x, y, z, x*y,
		x*z, y*z, x*x, y*y, z*z,
		x*y*z, x*x*x, x*y*y, x*z*z, x*x*y,
		y*y*y, y*z*z, x*x*z, y*y*z, z*z*z};
	for (int i = 0; i < 20; i++)
		m[i] = t[i];
}

// evaluate a polynomial of degree 3 in 3 variables
static long double eval_pol20l(long double c[20], long double X[3])
{
	long double m[20];
	eval_mon20(m, X);
	return scalar_product(m, c, 20);
}

// evaluate a rational function p(x)/q(x)
long double rpceval33(long double p[20], long double q[20], long double x[3])
{
	return eval_pol20l(p, x) / eval_pol20l(q, x);
}

// evaluate the error of a list of points by the given p,q
long double rpc33_error(long double p[20], long double q[20],
		long double (*x)[3], long double *f, int n)
{
	long double v[n];
	for (int i = 0; i < n; i++)
		v[i] = rpceval33(p, q, x[i]);

	long double r = 0;
	for (int i = 0; i < n; i++)
		r = pow(f[i] - v[i], 2);
	return 100 * sqrt(r/n);
}

// evaluate the error of a list of points by the given pq (sets q[0]=1)
long double rpc33_errorpq(long double pq[39],
		long double (*x)[3], long double *f, int n)
{
	long double q[20] = {1};
	for (int i = 1; i < 20; i++)
		q[i] = pq[19+i];
	return rpc33_error(pq, q, x, f, n);
}

// least-squares fit a 3rd degree RPC model to the given data values in 3D
long double rpcfit33(
		long double p[20],    // output p coefficients
		long double q[20],    // output q coefficients (q[0] = 1)
		long double (*x)[3],  // input data positions
		long double *f,       // input data values
		int n                 // number of input points
		)
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
	q[0] = 1;
	for (int i = 1; i < 20; i++)
		q[i] = 0;

	// build matrix "M" of evaluated monomials
	long double M[n][39];
	for (int i = 0; i < n; i++)
	{
		eval_mon20(M[i], x[i]);
		for (int j = 1; j < 20; j++)
			M[i][19+j] = -f[i] * M[i][j];
	}

	long double best_error_so_far = INFINITY;
	int number_of_iterations = 13;
	for (int iteration = 0; iteration < number_of_iterations; iteration++)
	{
		// vector "w" of weights
		long double w[n];
		for (int i = 0; i < n; i++)
			w[i] = pow( eval_pol20l(q, x[i]), -2);

		// matrix "A" of the system
		long double A[39][39] = {{0}};
		for (int i = 0; i < 39; i++)
		for (int j = 0; j < 39; j++)
		for (int k = 0; k < n ; k++)
			A[i][j] += w[k] * M[k][i] * M[k][j];

		// vector "b"
		long double b[39] = {0};
		for (int i = 0; i < 39; i++)
		for (int j = 0; j < n ; j++)
			b[i] += w[j] * f[j] * M[j][i];

		// solve the linear problem
		long double pq[39];
		solve_spd(pq, (void*)A, b, 39);

		// update if it improves error
		long double e = rpc33_errorpq(pq, x, f, n);
		if (e < best_error_so_far)
		{
			best_error_so_far = e;
			for (int i = 0; i < 20; i++)
				p[i] = pq[i];
			for (int i = 1; i < 20; i++)
				q[i] = pq[19+i];
		}
	}

	return best_error_so_far;
}
