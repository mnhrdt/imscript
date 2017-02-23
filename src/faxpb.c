#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr, "usage:\n\t%s a b [in [out]]\n", *v);
		//                          0 1 2  3   4
		return EXIT_FAILURE;
	}
	float a = atof(v[1]);
	float b = atof(v[2]);
	char *in = c > 3 ? v[3] : "-";
	char *out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	for (int i = 0; i < w * h * pd; i++)
		x[i] = a*x[i] + b;
	iio_write_image_float_vec(out, x, w, h, pd);
	free(x);
	return EXIT_SUCCESS;
}
