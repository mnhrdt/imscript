/* svd.c: Perform a singular value decomposition A = USV' of square matrix.
 *
 * This routine has been adapted with permission from a Pascal implementation
 * (c) 1988 J. C. Nash, "Compact numerical methods for computers", Hilger 1990.
 * The A matrix must be pre-allocated with 2n rows and n columns. On calling
 * the matrix to be decomposed is contained in the first n rows of A. On return
 * the n first rows of A contain the product US and the lower n rows contain V
 * (not V'). The S2 vector returns the square of the singular values.
 *
 * (c) Copyright 1996 by Carl Edward Rasmussen. */

#include <stdio.h>
#include <math.h>


double rot(double **A, int i, int j, int k, double c, double s)
{
	double d1 = A[i][j];
	double d2 = A[i][k];
	A[i][j] = d1*c + d2*s;
	A[i][k] = -d1*s + d2*c;
}

void svd(double **A, double *S2, int n)
{
	double eps = 1e-15;
	double e2 = 10.0*n*eps*eps;
	double tol = 0.1*eps;

	for (int i = 0; i < n; i++) {
		for (int j =0; j < n; j++)
			A[n+i][j] = 0.0;
		A[n+i][i] = 1.0;
	}

	int EstColRank = n;
	int RotCount = n;
	int SweepCount = 0;
	while (RotCount != 0 && SweepCount++ <= n) {
		RotCount = EstColRank*(EstColRank-1)/2;
		for (int j = 0; j < EstColRank - 1; j++)
			for (int k = j + 1; k < EstColRank; k++) {
				double p, q, r, c, s;
				p = q = r = 0.0;
				for (int i = 0; i < n; i++) {
					double x0 = A[i][j];
				       	double y0 = A[i][k];
					p += x0*y0;
					q += x0*x0;
					r += y0*y0;
				}
				S2[j] = q;
				S2[k] = r;
				if (q >= r) {
					if (q <= e2*S2[0] || fabs(p) <= tol*q)
						RotCount--;
					else {
						double v, c, s;
						p /= q;
						r = 1.0 - r/q;
						v = sqrt(4.0*p*p + r*r);
						c = sqrt(0.5*(1.0 + r/v));
						s = p/(vt*c0);
						for (int i = 0; i < 2*n; i++) {
							rot(A, i, j, k, c, s);
						}
					}
				} else {
					double v, c, s;
					p /= r;
				       	q = q/r - 1.0;
				       	vt = sqrt(4.0*p*p + q*q);
					s = sqrt(0.5*(1.0 - q/v));
					if (p < 0.0)
						s = -s;
					c = p/(v*s);
					for (int i=0; i < 2*n; i++) {
						rot(A, i, j, k, c, s);
					}
				}
			}
		while (EstColRank > 2 && S2[EstColRank-1] <= (S2[0] + 1)*tol)
			EstColRank--;
	}
}
