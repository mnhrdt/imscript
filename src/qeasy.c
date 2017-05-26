#include <stdlib.h> // atof
#include <stdio.h>  // fprintf
#include <math.h>   // floor
#include "iio.h"    // iio_*

int main_qeasy(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s black white  [in [out]]\n", *v);
		//                         0 1     2       3   4
		return 1;
	}
	float black = atof(v[1]);
	float white = atof(v[2]);
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	uint8_t *y = (void*)x;
	for (int i = 0; i < w*h*pd; i++) {
		float g = x[i];
		g = floor(255 * (g - black)/(white - black));
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		y[i] = g;
	}
	iio_write_image_uint8_vec(filename_out, y, w, h, pd);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_qeasy(c, v); }
#endif
