
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
static void solve_spdl(long double *x, long double *A, long double *b, int n)
{
	// factor the matrix A = L * L'
	long double L[n*n];
	cholesky(L, A, n);

	// solve the two triangular systems L * L' * x = b
	long double t[n];
	solve_lower(t, L, b, n);
	solve_upper(x, L, t, n);
}

// solve a symmetric, positive definite system
static void solve_spdf(float *xf, float *Af, float *bf, int n)
{
	long double x[n], A[n*n], b[n];
	for (int i = 0; i < n; i++)
		b[i] = bf[i];
	for (int i = 0; i < n*n; i++)
		A[i] = Af[i];
	solve_spdl(x, A, b, n);
	for (int i = 0; i < n; i++)
		xf[i] = x[i];
}

// solve a symmetric, positive definite system
static void solve_spd(double *xf, double *Af, double *bf, int n)
{
	long double x[n], A[n*n], b[n];
	for (int i = 0; i < n; i++)
		b[i] = bf[i];
	for (int i = 0; i < n*n; i++)
		A[i] = Af[i];
	solve_spdl(x, A, b, n);
	for (int i = 0; i < n; i++)
		xf[i] = x[i];
}
