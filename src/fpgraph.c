#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

static double fpenergy(double p, double *x, int n, double m)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += pow(fabs(x[i] - m), p);
	return r/n;
}

static double fpenergy_normalized(double p, double *x, int n, double m)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += pow(fabs(x[i] - m), p);
	return pow(r/n,1/p);
}

#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	bool normalize = pick_option(&c, &v, "n", 0);
	if (c <= 6)
		return fprintf(stderr, "usage:\n\t"
				"%s p xmin xmax nsamples x1 ... xN\n", *v);
	//                        0 1 2    3    4        N+5
	int n = c - 5;
	double x[n];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[5+i]);
	double p = atof(v[1]);
	double xmin = atof(v[2]);
	double xmax = atof(v[3]);
	int nsamples = atoi(v[4]);

	for (int i = 0; i < nsamples; i++)
	{
		double m = xmin + i * (xmax - xmin) / (nsamples-1);
		double e = fpenergy(p, x, n, m);
		if (normalize)
			e = fpenergy_normalized(p, x, n, m);
		printf("%g %g\n", m, e);
	}

	return 0;
}
