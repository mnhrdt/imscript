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

/* This function has been modified by Enric Meinhardt-Llopis to accept an
 * alternative calling convention */

#include <stdio.h>
#include <math.h>


void svd(double *out_V, double *out_S2, double *in_A, int n)
{
	double (*A)[n] = (void*)in_A;
	double (*V)[n] = (void*)out_V;
	double *S2 = out_S2;
	int  i, j, k, EstColRank = n, RotCount = n, SweepCount = 0,
	     slimit = n;//(n<120) ? 30 : n/4;
	double eps = 1e-15, e2 = 10.0*n*eps*eps, tol = 0.1*eps, vt, p, x0,
	       y0, q, r, c0, s0, d1, d2;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			V[i][j] = 0.0;
		V[i][i] = 1.0;
	}

	while (RotCount != 0 && SweepCount++ <= slimit) {

		RotCount = EstColRank*(EstColRank - 1)/2;

		for (j = 0; j < EstColRank - 1; j++)
			for (k = j+1; k < EstColRank; k++) {
				p = q = r = 0.0;
				for (i = 0; i < n; i++) {
					x0 = A[i][j];
					y0 = A[i][k];
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
						p /= q;
						r = 1.0 - r/q;
						vt = sqrt(4.0*p*p + r*r);
						c0 = sqrt(0.5*(1.0 + r/vt));
						s0 = p/(vt*c0);
						for (i = 0; i < n; i++) {
							d1 = A[i][j];
							d2 = A[i][k];
							A[i][j] = d1*c0+d2*s0;
							A[i][k] = -d1*s0+d2*c0;
							d1 = V[i][j];
							d2 = V[i][k];
							V[i][j] = d1*c0+d2*s0;
							V[i][k] = -d1*s0+d2*c0;
						}
					}
				} else {
					p /= r;
					q = q/r - 1.0;
					vt = sqrt(4.0*p*p + q*q);
					s0 = sqrt(0.5*(1.0 - q/vt));
					if (p < 0.0)
						s0 = -s0;
					c0 = p/(vt*s0);
					for (i = 0; i < n; i++) {
						d1 = A[i][j];
						d2 = A[i][k];
						A[i][j] = d1*c0+d2*s0;
						A[i][k] = -d1*s0+d2*c0;
						d1 = V[i][j];
						d2 = V[i][k];
						V[i][j] = d1*c0+d2*s0;
						V[i][k] = -d1*s0+d2*c0;
					}
				}
			}
		while (EstColRank > 2 && S2[EstColRank-1] <= (S2[0] + tol)*tol)
			EstColRank--;
	}
	//if (SweepCount > slimit)
	//	fprintf(stderr, "Warning: "
	//			"Reached maximum number of sweeps (%d) "
	//			"in SVD routine...\n" ,slimit);
}
