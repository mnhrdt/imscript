#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	if (pd == 1 || pd == 3)
		iio_save_image_float_vec(out, x, w, h, pd);
	else if (pd == 2) {
		for (int i = 0; i < w*h; i++)
			x[i] = x[2*i];
		iio_save_image_float_vec(out, x, w, h, 1);
	} else if (pd == 4) {
		for (int i = 0; i < w*h; i++) {
			x[3*i] = x[4*i];
			x[3*i+1] = x[4*i+1];
			x[3*i+2] = x[4*i+2];
		}
		iio_save_image_float_vec(out, x, w, h, 3);
	} else
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
