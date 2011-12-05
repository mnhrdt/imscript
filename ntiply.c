#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 2 && c != 4 && c != 3) {
		fprintf(stderr, "usage:\n\t%s n [in [out]]\n", *v);
		//                          0 1  2   3
		return EXIT_FAILURE;
	}
	int n = atof(v[1]);
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	float *y = malloc(n*n*w*h*pd*sizeof*y);
	for (int j = 0; j < n*h; j++)
	for (int i = 0; i < n*w; i++)
	for (int l = 0; l < pd; l++)
		y[pd*(n*w*j+i)+l] = x[pd*(w*(j/n)+i/n)+l];
	iio_save_image_float_vec(out, y, n*w, n*h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
