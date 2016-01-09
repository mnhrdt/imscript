// downwards neighbor interpolation

#include <math.h>

void unshadow(float *xx, int w, int h)
{
	float (*x)[w] = (void*)xx;
	for (int i = 0; i < w; i++)
		for (int j = h - 2; j >= 0; j--)
			if (isnan(x[j][i]) && !isnan(x[j+1][i]))
				x[j][i] = x[j+1][i];
}



#ifndef OMIT_UNSHADOW_MAIN
#define USE_UNSHADOW_MAIN
#endif

#ifdef USE_UNSHADOW_MAIN
#include <stdio.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *filename_mask = pick_option(&c, &v, "m", "");
	int help_argument = (int)pick_option(&c, &v, "h", 0);
	if (help_argument || (c != 1 && c != 2 && c != 3)) {
		fprintf(stderr, "usage:\n\t%s [in.tiff [out.tiff]]\n", *v);
		//                          0  1        2
		return 1;
	}
	char *filename_in   = c > 1 ? v[1] : "-";
	char *filename_out  = c > 2 ? v[2] : "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	if (filename_mask && *filename_mask) {
		int mw, mh;
		float *m = iio_read_image_float(filename_mask, &mw, &mh);
		for (int i = 0; i < mw*mh; i++)
			if (i < w*h && m[i])
				x[i] = NAN;
	}


	unshadow(x, w, h);

	iio_save_image_float(filename_out, x, w, h);

	return 0;
}
#endif//USE_UNSHADOW_MAIN
