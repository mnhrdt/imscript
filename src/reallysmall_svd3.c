// ISO C function for computing A=USV'
//
// Based on the 1988 book by J.C.Nash "Compact numerical methods for computers"
// Based on code found online by Carl Edward Rasmussen
//
// INPUT
// in_A: a square positive semi-definite symmetric matrix
// n: dimension of the input matrix
//
// OUTPUT:
// out_V: filled with the matrix V
// out_S2: filled with the squares of the singular values
// in_A: overwritten with the product US
//
#include <math.h>
void svd(double *out_V, double *out_S2, double *in_A, int n)
{
	double (*A)[n] = (void*)in_A;
	double (*V)[n] = (void*)out_V;
	double *S2 = out_S2;

	double eps = 1e-15;
	double e2 = 10.0*n*eps*eps;
	double tol = 0.1*eps;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			V[i][j] = 0.0;
		V[i][i] = 1.0;
	}

	int EstColRank = n;
	int RotCount = n;
	int SweepCount = 0;

	while (RotCount != 0 && SweepCount++ <= n) {

		RotCount = EstColRank*(EstColRank - 1)/2;

		for (int j = 0; j < EstColRank - 1; j++)
			for (int k = j + 1; k < EstColRank; k++) {
				double p, q, r;
				p = q = r = 0.0;
				for (int i = 0; i < n; i++) {
					double x = A[i][j];
					double y = A[i][k];
					p += x*y;
					q += x*x;
					r += y*y;
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
						s = p/(v*c);
						for (int i = 0; i < n; i++) {
							double a, b;
							a = A[i][j];
							b = A[i][k];
							A[i][j] = a*c + b*s;
							A[i][k] = -a*s + b*c;
							a = V[i][j];
							b = V[i][k];
							V[i][j] = a*c + b*s;
							V[i][k] = -a*s + b*c;
						}
					}
				} else {
					double v, c, s;
					p /= r;
					q = q/r - 1.0;
					v = sqrt(4.0*p*p + q*q);
					s = sqrt(0.5*(1.0 - q/v));
					if (p < 0.0)
						s = -s;
					c = p/(v*s);
					for (int i = 0; i < n; i++) {
						double a, b;
						a = A[i][j];
						b = A[i][k];
						A[i][j] = a*c + b*s;
						A[i][k] = -a*s + b*c;
						a = V[i][j];
						b = V[i][k];
						V[i][j] = a*c + b*s;
						V[i][k] = -a*s + b*c;
					}
				}
			}
		while (EstColRank > 2 && S2[EstColRank-1] <= (S2[0] + tol)*tol)
			EstColRank--;
	}
}
