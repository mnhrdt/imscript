#include <math.h> // sin, cos
#include <stdio.h> // fprintf

#include "random.c"

static void radar_sim_horizontal(
		float *y,      // output: radar simulation
		float *x,      // input: height raster in pixel units
		int w,         // raster width
		int h,         // raster height
		float a        // angle pi=left pi/2=top 0=right
		)
{
	// TODO: remove invisible parts (preprocess with shadowcast)
	// TODO: find reflectors and simulate them
	//
	fprintf(stderr, "radarsim a=%g\n", a);
	for (int i = 0; i < w*h; i++)
		y[i] = 100;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p = w/2 + cos(a) * (i-w/2) - sin(a) * x[j*w+i];
		int ip = p;
		float fp = p - ip;
		//float r = 1; // may be random number
		//float r = random_uniform();
		//float r0 = random_uniform();
		if (ip >= 0 && ip < w+1)
		{
			y[j*w+ip+0] += (1-fp)*random_uniform();
			y[j*w+ip+1] += fp*random_uniform();
		}
		else
			if (j == 10)
				fprintf(stderr, "bad ip=%d\n", ip);
	}
}

#ifndef OMIT_MAIN_SARSIM
#define MAIN_SARSIM
#endif


#ifdef MAIN_SARSIM
#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.c"
#include "iio.h"      // library for image input/output
int main(int c, char *v[])
{
	// process input arguments
	if (c < 2 || c > 4) {
		fprintf(stderr, "usage:\n\t%s a [dem_in [dem_out]]\n", *v);
		//                          0 1  2       3
		return 1;
	}
	float param_a = atof(v[1]); // angle
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// read input dem
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *y = xmalloc(w * h * sizeof*x);

	// simulate radar
	radar_sim_horizontal(y, x, w, h, param_a*M_PI/180);

	// save the output image
	iio_write_image_float(filename_out, y, w, h);

	// cleanup (unnecessary) and exit
	free(y);
	free(x);
	return 0;
}
#endif//MAIN_SARSIM
