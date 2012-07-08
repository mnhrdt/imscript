#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "parsenumbers.c"

#include "vvector.h"
static void invert_homography(double invH[9], double H[9])
{
	double h[3][3] = { {H[0], H[1], H[2]},
			{H[3], H[4], H[5]},
			{H[6], H[7], H[8]}};
	double det, ih[3][3];
	INVERT_3X3(ih, det, h);
	for (int i = 0; i < 9; i++) invH[i] = ih[0][i];
}


static void projective_map(double y[2], double H[9], double x[2])
{
	double z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
}

#include "smapa.h"
SMART_PARAMETER(HOMI,0)

int main(int c, char *v[])
{
	if (c != 10) {
		fprintf(stderr, "usage:\n\t%s h1 ... h9 < points \n", *v);
		//                          0 1      9
		return EXIT_FAILURE;
	}
	double H[9];
	for (int i = 0; i < 9; i++)
		H[i] = atof(v[1+i]);
	if (HOMI() > 0)
		invert_homography(H, H);

	int n;
	double *t = read_ascii_doubles(stdin, &n);
	n /= 2;
	for (int i = 0; i < n; i++)
		projective_map(t+2*i, H, t+2*i);
	for (int i = 0; i < n; i++)
		printf("%g %g\n", t[2*i], t[2*i+1]);
	return EXIT_SUCCESS;
}
