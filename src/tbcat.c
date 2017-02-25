#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "smapa.h"
#include "xmalloc.c"
#include "getpixel.c"
#include "pickopt.c"

#define BAD_MAX(a,b) (a)<(b)?(b):(a);

SMART_PARAMETER(BACKGROUND,0)

int main_tbcat(int c, char *v[])
{
	if (c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s a b [ab]\n", *v);
		//                          0 1 2  3
		return EXIT_FAILURE;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *filename_out = c > 3 ? v[3] : "-";

	int w[3], h[3], pd[3];
	float *x = iio_read_image_float_vec(filename_a, w+0, h+0, pd+0);
	float *y = iio_read_image_float_vec(filename_b, w+1, h+1, pd+1);
	w[2] = BAD_MAX(w[0],w[1]);
	h[2] = h[0] + h[1];
	pd[2] = BAD_MAX(pd[0], pd[1]);

	float *z = xmalloc(w[2] * h[2] * pd[2] * sizeof*y);
	for (int i = 0; i < w[2] * h[2] * pd[2]; i++)
		z[i] = BACKGROUND();

	for (int j = 0; j < h[0]; j++)
	for (int i = 0; i < w[0]; i++)
		for (int l = 0; l < pd[2]; l++) {
			float s = getsample_1(x, w[0], h[0], pd[0], i, j, l);
			setsample_0(z, w[2], h[2], pd[2], i, j, l, s);
		}

	for (int j = 0; j < h[1]; j++)
	for (int i = 0; i < w[1]; i++)
		for (int l = 0; l < pd[2]; l++) {
			float s = getsample_1(y, w[1], h[1], pd[1], i, j, l);
			setsample_0(z, w[2], h[2], pd[2], i, j+h[0], l, s);
		}

	iio_write_image_float_vec(filename_out, z, w[2], h[2], pd[2]);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_tbcat(c, v); }
#endif
