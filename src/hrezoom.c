#include <stdio.h>
#include <stdlib.h>

#include "vvector.h"

int main(int c, char *v[])
{
	if (c != 11) {
		fprintf(stderr, "usage:\n\t%s a b p c d q r s t zoom\n", *v);
		//                          0 1 2 3 4 5 6 7 8 9 10
		return 1;
	}
	double a[3][3], z = atof(v[10]);
	for (int i = 0; i < 9; i++)
		a[0][i] = atof(v[1+i]);
	double zp[3][3] = {{1/z, 0, 0}, {0, 1/z, 0}, {0,0,1}};
	double zm[3][3] = {{z, 0, 0}, {0, z, 0}, {0,0,1}};
	double oa[3][3];
	MATRIX_PRODUCT_3X3(oa,a,zp);
	MATRIX_PRODUCT_3X3(a,zm,oa);
	for (int i = 0; i < 9; i++)
		printf(" %lf", a[0][i]);
	printf("\n");
	return 0;
}
