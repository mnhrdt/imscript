#include <stdio.h>
#include <stdlib.h>
#include "iio.h"


void projective_map(double y[2], double A[9], double x[2])
{
	double p = A[0]*x[0] + A[1]*x[1] + A[2];
	double q = A[3]*x[0] + A[4]*x[1] + A[5];
	double y2 = A[6]*x[0] + A[7]*x[1] + A[8];
	y[0] = p/y2;
	y[1] = q/y2;
}

#include "xmalloc.c"
#include "parsenumbers.c"
#include "pickopt.c"

int main(int c, char **v)
{
	char *hlstring = pick_option(&c, &v, "l", "1 0 0 0 1 0 0 0 1");
	char *hrstring = pick_option(&c, &v, "r", "1 0 0 0 1 0 0 0 1");
	if (c != 1) {
		fprintf(stderr, "usage:\n\t%s [-l \"h1...h9\"] [-r \"h1...h9\"] <pairs.txt >hpairs.txt\n", *v);
		return 0;
	}

	double HA[9], HB[9];
	int na = parse_doubles(HA, 9, hlstring);
	int nb = parse_doubles(HB, 9, hrstring);
	if (na != 9 || nb != 9)
		fail("%d,%d != 9,9", na, nb);

	int n;
	double *p = read_ascii_doubles(stdin, &n);
	n /= 4;

	for (int i = 0; i < n; i++)
	{
		double *x = p + 4*i;
		double *y = p + 4*i + 2;
		projective_map(x, HA, x);
		projective_map(y, HB, y);

		printf("%g %g %g %g\n", x[0], x[1], x[2], x[3]);
	}



	free(p);
	return 0;
}
