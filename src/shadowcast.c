#include <math.h>
#include <stdio.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void cast_vertical_shadows(float *xx, int w, int h, float alpha)
{
	// pointer fo easy access ( x[j][i] == xx[j*w+i] )
	float (*x)[w] = (void*)xx;

	// compute slope of the rays
	float slope = tan(alpha * M_PI / 180);
	fprintf(stderr, "casting shadows with slope %g\n", slope);

	// process each column independently
	for (int i = 0; i < w; i++)
	{
		int l = 0;
		for (int j = 0; j < h; j++)
			if (x[j][i] < slope * (j - l) + x[l][i])
				x[j][i] = NAN;
			else
				l = j;
	}
}

#define MAIN_VERTSHADOW
#ifdef MAIN_VERTSHADOW
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"      // library for image input/output
int main(int c, char *v[])
{
	// process input arguments
	if (c < 2 || c > 4) {
		fprintf(stderr, "usage:\n\t%s alpha [dem_in [dem_out]]\n", *v);
		//                          0 1      2       3
		return 1;
	}
	float alpha = atof(v[1]);
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (pd != 1)
		return fprintf(stderr, "1D-valud image expected\n");

	// cast the vertical shadows
	cast_vertical_shadows(x, w, h, alpha);

	// save the output image
	iio_save_image_float_vec(filename_out, x, w, h, 1);

	// cleanup (unnecessary) and exit
	return 0;
}
#endif//MAIN_VERTSHADOW

