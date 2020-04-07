#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"

int main(int c, char *v[])
{
	_Bool t = pick_option(&c, &v, "t", 0);
	_Bool T = pick_option(&c, &v, "T", 0);
	if (c != 3)
		return fprintf(stderr,
				"usage:\n\t%s [-t] in.tif out_%%d.tif\n", *v);
		//                          0       1      2
	char *filename_in = v[1];
	char *filepat_out = v[2];
	int w, h, pd;
	if (!t && !T) {
		fprintf(stderr, "getbands standard...\n");
		float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
		fprintf(stderr, "getbands standard w=%d h=%d pd=%d\n",w,h,pd);
		for (int i = 0; i < pd; i++)
		{
			char n[FILENAME_MAX];
			snprintf(n, FILENAME_MAX, filepat_out, i);
			iio_write_image_float(n, x + w*h*i, w, h);
		}
	} else if (!T) {
		float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
		for (int i = 0; i < h; i++)
		{
			char n[FILENAME_MAX];
			snprintf(n, FILENAME_MAX, filepat_out, i);
			iio_write_image_float(n, x + w*pd*i, pd, w);
		}
	} else {
		float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
		float *y = malloc(h*pd*sizeof*y);
		for (int i = 0; i < h*pd; i++)
			y[i] = 42;
		fprintf(stderr, "getbands T option w=%d h=%d pd=%d\n",w,h,pd);
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			for (int l = 0; l <pd; l++)
				y[j*pd+l] = x[(j*w+i)*pd+l];
			char n[FILENAME_MAX];
			snprintf(n, FILENAME_MAX, filepat_out, i);
			iio_write_image_float(n, y, pd, h);
		}
	}
	return 0;
}
