#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "xmalloc.c"
#include "getpixel.c"
#include "pickopt.c"

#define BAD_MAX(a,b) (a)<(b)?(b):(a);

#include "smapa.h"
SMART_PARAMETER_SILENT(BACKGROUND,0)

int main_lrcat_two(int c, char *v[])
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
	w[2] = w[0] + w[1];
	h[2] = BAD_MAX(h[0], h[1]);
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
			setsample_0(z, w[2], h[2], pd[2], i+w[0], j, l, s);
		}

	iio_write_image_float_vec(filename_out, z, w[2], h[2], pd[2]);
	return EXIT_SUCCESS;
}

int main_lrcat(int c, char *v[])
{
	char *filename_out = pick_option(&c, &v, "o", "-");
	int n = c - 1;
	char *filename[n+1];
	for (int i = 0; i < n; i++)
		filename[i] = v[i+1];
	filename[n] = filename_out;

	int w[n+1], h[n+1], pd[n+1];
	float *x[n+1];
	w[n] = h[n] = pd[n] = 0;
	for (int k = 0; k < n; k++)
	{
		x[k] = iio_read_image_float_vec(filename[k], w+k, h+k, pd+k);
		w[n] += w[k];
		h[n] = BAD_MAX(h[n], h[k]);
		pd[n] = BAD_MAX(pd[n], pd[k]);
	}
	x[n] = xmalloc(w[n] * h[n] * pd[n] * sizeof(x[0][0]));
	for (int i = 0; i < w[n] * h[n] * pd[n]; i++)
		x[n][i] = BACKGROUND();

	int ok = 0;
	for (int k = 0; k < n; k++)
	{
		for (int j = 0; j < h[k]; j++)
		for (int i = 0; i < w[k]; i++)
		for (int l = 0; l < pd[2]; l++) {
			float s = getsample_1(x[k], w[k], h[k], pd[k], i, j, l);
			setsample_0(x[n], w[n], h[n], pd[n], i+ok, j, l, s);
		}
		ok += w[k];
	}

	iio_write_image_float_vec(filename[n], x[n], w[n], h[n], pd[n]);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_lrcat(c, v); }
#endif
