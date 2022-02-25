// slice along t "smile" image

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 4)
		return fprintf(stderr,
				"usage:\n\t%s guide in out\n", *v);
		//                          0 1     2  3
	char *filename_guide = v[1];
	char *filename_in    = v[2];
	char *filename_out   = v[3];
	int w, h, pd, ww, hh;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *g = iio_read_image_float(filename_guide, &ww, &hh);
	if (ww != w || hh != h)
		exit(fprintf(stderr, "guide/input mismatch %dx%d != %dx%d\n",
					ww, hh, w, h));
	float *y = malloc(w*h*sizeof*y);
	for (int i = 0; i < w*h; i++)
	{
		int j = g[i];
		if (j < 0 || j >= pd)
			y[i] = NAN;
		else
			y[i] = x[i*pd+j];
	}
	iio_write_image_float_vec(filename_out, y, w, h, 1);
	return 0;
}
