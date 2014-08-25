#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "smapa.h"
#include "random.c"
#include "iio.h"

SMART_PARAMETER_SILENT(SRAND, 0)

int main(int c, char *v[])
{
	if (c != 4 && c != 3 && c != 2) {
		fprintf(stderr, "usage:\n\t%s sigma [in [out]]\n", *v);
		//                          0 1      2   3
		return EXIT_FAILURE;
	}
	float s = atof(v[1]);
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);

	xsrand(SRAND());
	for (int i = 0; i < w*h*pd; i++)
		x[i] += s * random_normal();
	iio_save_image_float_vec(out, x, w, h, pd);

	return EXIT_SUCCESS;
}
