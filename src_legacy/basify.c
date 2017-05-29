// affine projection matrix determined by a satellite direction

#include <assert.h>
#include <math.h>

static double norm(double a, double b, double c)
{
	return hypot(hypot(a, b), c);
}

void basify(double out_P[8], double in_p[3])
{
	// P = normalized input vector
	double p = in_p[0];
	double q = in_p[1];
	double r = in_p[1];
	double np = norm(p, q, r);
	p /= np; q /= np; r /= np;

	// A = orthogonal vector to P, to be computed
	double a, b, c;

	double e = norm(p, q, r - 1);
	if (e < 1e-5) { // P is nearly vertical
		        //(needed to avoid gimbal lock numerics)
		a = r;
		b = 0;
		c = -q;
	} else { // general case
		a = -q;
		b = p;
		c = 0;
	}

	// normalize vector A
	double na = norm(a, b, c);
	a /= na; b /= na; c /= na;

	// X = P x A
	double x = q * c - r * b;
	double y = r * a - p * c;
	double z = p * b - q * a;

	// normalize vector X (unnneeded)
	double nx = norm(x, y, z);
	assert(fabs(nx - 1) < 1e-5);
	//x /= nx; y /= nx; z /= nx;
	
	// 


}


#include <stdio.h>
#include <stdlib.h>
int main(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s p q r > P.txt\n");
		return 1;
	}
	double p[3] = { atof(v[1]), atof(v[2]), atof(v[3]) };
	double P[8];

	basify(P, p);

	for (int i = 0; i < 8; i++)
		printf("%g%c", P[i], i==7?'\n':' ');
	
	return 0;
}
