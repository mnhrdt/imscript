#include <stdio.h>
#include <stdlib.h>
#include "iio.h"


int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

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

	qsort(x, w*h*pd, sizeof*x, compare_floats);

	iio_write_image_float_vec(out, x, w, h, pd);
	free(x);
	return EXIT_SUCCESS;
}
