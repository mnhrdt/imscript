#include <math.h>

static double dot(double *x, double *y, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

static double norm(double *x, double n)
{
	return sqrt(dot(x, x, n));
}

void gram_schmidt(double *yy, double *xx, int n, int m)
{
	double (*x)[n] = (void*)xx;
	double (*y)[n] = (void*)yy;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
			y[j][i] = x[j][i];
		for (int k = 0; k < j; k++)
		{
			double yjyk = dot(y[j], y[k], n);
			double ykyk = dot(y[k], y[k], n);
			for (int i = 0; i < n; i++)
				y[j][i] -= (yjyk / ykyk) * y[k][i];
		}
		double nxj = norm(x[j], n);
		for (int i = 0; i < n; i++)
			y[j][i] = x[j][i] / nxj;
	}
}

//// gram schmidt orthogonalization of m vectors in R^n
//void gram_schmidt_verbose(
//		double *yy,  // output (m orthogonal vectors of length n)
//		double *xx,  // input  (m independent vectors of length n)
//		int n,       // dimension
//		int m        // number of input vectors
//		)
//{
//	// pointers and macros for easy algebra
//	double (*x)[n] = (void*)xx;     // xx[j*n+i] == x[j][i]
//	double (*y)[n] = (void*)yy;     // yy[j*n+i] == x[j][i]
//#define FORALL for (int i=0; i<n; i++)  // for symbolic index notation
//
//	for (int j = 0; j < m; j++) // j = number of basis vectors already built
//	{
//		// start by copying the next vector
//		FORALL  y[j][i] = x[j][i];
//
//		// then substract from it each projection
//		for (int k = 0; k < j; k++)
//		{
//			double yjxk = dot(y[j], x[k], n);
//			double xkxk = dot(x[k], x[k], n);
//			FORALL  y[j][i] -= (yjxk / xkxk) * x[k][i];
//		}
//
//		// finally normalize the result
//		FORALL  y[j][i] = x[j][i] / norm(x[j], n);
//	}
//}
//
