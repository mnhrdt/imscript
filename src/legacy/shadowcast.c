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

	if (alpha <= 0)
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
	else
		for (int i = 0; i < w; i++)
		{
			float slop = tan((180-alpha)*M_PI/180);
			int l = h-1;
			for (int j = h-1; j >= 0; j--)
				if (x[j][i] < slop * (l - j) + x[l][i])
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
#include "pickopt.c"  // function "pick_option" for processing args
int main(int c, char *v[])
{
	// process input arguments
	_Bool m = pick_option(&c, &v, "m", NULL);
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
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	if (pd != 1) {
		for (int i = 0; i < w*h; i++)
			for (int l = 1; l < pd; l++)
				x[i] += x[i+w*h*l];
		for (int i = 0; i < w*h; i++)
			x[i] /= pd;
	}

	// cast the vertical shadows
	cast_vertical_shadows(x, w, h, alpha);

	// if mask is requested, create a binary mask
	if (m) for (int i = 0; i < w*h; i++)
		x[i] = 255*isnan(x[i]);

	// save the output image
	iio_write_image_float_split(filename_out, x, w, h, 1);

	// cleanup (unnecessary) and exit
	return 0;
}
#endif//MAIN_VERTSHADOW

