#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 6 && c != 5 && c != 3) {
		fprintf(stderr, "usage:\n\t%s a b in1 in2 [out]\n", *v);
		//                          0 1 2 3   4    5
		return EXIT_FAILURE;
	}
	float a = atof(v[1]);
	float b = atof(v[2]);
	char *outfile = c > 5 ? v[5] : "-";

	int w[2], h[2], pd[2];
	float *x1 = iio_read_image_float_vec(v[3], w+0, h+0, pd+0);
	float *x2 = iio_read_image_float_vec(v[4], w+1, h+1, pd+1);
	if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1]) {
		fprintf(stderr, "images size mismatch\n");
		return EXIT_FAILURE;
	}
	for (int i = 0; i < *w * *h * *pd; i++)
		x1[i] = a*x1[i] + b*x2[i];
	iio_save_image_float_vec(outfile, x1, *w, *h, *pd);
	return EXIT_SUCCESS;
}
