#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s in.tif out_%%d.tif\n", *v);
		//                                 0 1      2
	char *filename_in = v[1];
	char *filepat_out = v[2];
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	for (int i = 0; i < pd; i++)
	{
		char n[FILENAME_MAX];
		snprintf(n, FILENAME_MAX, filepat_out, i);
		iio_write_image_float(n, x + w*h*i, w, h);
	}
	return 0;
}
