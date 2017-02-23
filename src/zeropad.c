#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "xmalloc.c"

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s mx my [in [out]]\n", *v);
		//                          0 1  2  3   4
		return EXIT_FAILURE;
	}
	int mx = atoi(v[1]);
	int my = atoi(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x_raw = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *y_raw = xmalloc( (w + mx)*(h + my)*pd*sizeof*y_raw);
	float (*x)[w][pd] = (void*)x_raw;
	float (*y)[w+mx][pd] = (void*)y_raw;

	for (int j = 0; j < h+my; j++)
	for (int i = 0; i < w+mx; i++)
	for (int l = 0; l < pd; l++)
		y[j][i][l] = (i<w && j<h) ? x[j][i][l] : 0;

	iio_write_image_float_vec(filename_out, y_raw, w+mx, h+my, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
