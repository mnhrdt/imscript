#include <stdio.h>
#include <stdlib.h>

#include "svd.c"

static void printmatrix(double *x, int nr, int nc, char *name)
{
	if (name)
		printf("%s = \n", name);
	for (int j = 0; j < nr; j++)
		for (int i = 0; i < nc; i++)
			printf("%- 12.4g%c", x[j*nc+i], i==(nc-1)?'\n':' ');
	if (name)
		printf("\n");
}


int main(int c, char *v[])
{
	if (c != 10) {
		fprintf(stderr, "usage:\n\t%s f1 ... f9\n", *v);
		//                          0 1  ... 9
		return 0;
	}
	double fm[9];
	for (int i = 0; i < 9; i++)
		fm[i] = atof(v[1+i]);

	double A[9][12] = {{0}};
	for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 3; k++)
			A[3*j + k][3*i + k] = fm[3*j + i];
	A[1][9] = A[2][10] = A[5][11] = -1;
	A[3][9] = A[6][10] = A[7][11] = 1;

	printmatrix(A[0], 9, 12, "A");

	// n=9, m=12
	double d[9], V[12*12], U[9*9];
	int r = svd_double(d, A[0], V, 12, U, 9);

	printmatrix(d, 9, 1, "d");
	printmatrix(U, 9, 9, "U");
	printmatrix(V, 12, 12, "V");


	return 0;
}
