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

static void check_solution(double f[9], double x[12])
{
	double FH[9], Fp[9];
	FH[0] = f[0]*x[0] + f[1]*x[3] + f[2]*x[6];
	FH[1] = f[0]*x[1] + f[1]*x[4] + f[2]*x[7];
	FH[2] = f[0]*x[2] + f[1]*x[5] + f[2]*x[8];

	FH[3] = f[3]*x[0] + f[4]*x[3] + f[5]*x[6];
	FH[4] = f[3]*x[1] + f[4]*x[4] + f[5]*x[7];
	FH[5] = f[3]*x[2] + f[4]*x[5] + f[5]*x[8];

	FH[6] = f[6]*x[0] + f[7]*x[3] + f[8]*x[6];
	FH[7] = f[6]*x[1] + f[7]*x[4] + f[8]*x[7];
	FH[8] = f[6]*x[2] + f[7]*x[5] + f[8]*x[8];

	Fp[0] = Fp[4] = Fp[8] = 0;
	Fp[1] = x[9];
	Fp[3] = -x[9];
	Fp[2] = x[10];
	Fp[6] = -x[10];
	Fp[5] = x[11];
	Fp[7] = -x[11];

	printmatrix(FH, 3, 3, "FH");
	printmatrix(Fp, 3, 3, "Fp");
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

	// define the linear system
	double A[9][12] = {{0}};
	for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 3; k++)
			A[3*j + k][3*i + k] = fm[3*j + i];
	A[1][9] = A[2][10] = A[5][11] = -1;
	A[3][9] = A[6][10] = A[7][11] = 1;

	printmatrix(A[0], 9, 12, "A");

	// solve the system by SVD
	// n=9, m=12
	double d[9], V[12*12], U[9*9];
	int r = svd_double(d, A[0], V, 12, U, 9);
	printf("SVDr = %d\n", r);
	if (r) abort();

	printmatrix(d, 1, 9, "d");
	printmatrix(U, 9, 9, "U");
	printmatrix(V, 12, 12, "V");

	double x[4][12];
	for (int j = 0; j < 4; j++)
	{
		double e[12] = {0};
		int k = 12 - j - 1;
		e[k] = 1;
		for (int i = 0; i < 12; i++)
			x[j][i] = V[12*i+k];
		printf("solution number %d\n", j);
		printmatrix(x[j], 1, 12, "xj");
	}

	for (int j = 0; j < 4; j++)
	{
		printf("check solution number %d\n", j);
		check_solution(fm, x[j]);
	}



	return 0;
}
