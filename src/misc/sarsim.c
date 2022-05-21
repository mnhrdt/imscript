#include <assert.h>
#include <math.h> // sin, cos, hypot
#include <stdio.h> // fprintf
#include <string.h> // memcpy

#include "random.c"

static void cast_horizontal_shadows(float *xx, int w, int h, float alpha)
{
	// pointer fo easy access ( x[j][i] == xx[j*w+i] )
	float (*x)[w] = (void*)xx;

	// compute slope of the rays
	float slope = tan(alpha * M_PI / 180);
	//fprintf(stderr, "casting shadows with slope %g\n", slope);

	if (alpha <= 0)
		// process each row independently
		for (int j = 0; j < h; j++)
		{
			int l = 0;
			for (int i = 0; i < w; i++)
				if (x[j][i] < slope * (i - l) + x[j][l])
					x[j][i] = NAN;
				else
					l = i;
		}
	else
		for (int j = 0; j < h; j++)
		{
			float slop = tan((180-alpha)*M_PI/180);
			int l = w-1;
			for (int i = w-1; i >= 0; i--)
				if (x[j][i] < slop * (l - i) + x[j][l])
					x[j][i] = NAN;
				else
					l = i;
		}
}

static void radar_sim_horizontal(
		float *y,      // output: radar simulation
		float *x,      // input: height raster in pixel units
		int w,         // raster width
		int h,         // raster height
		float a        // angle pi=left pi/2=top 0=right
		)
{
	fprintf(stderr, "radarsim a=%g\n", a);
	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	// TODO: remove invisible parts (preprocess with shadowcast)
	float *t = malloc(w*h*sizeof*t);
	memcpy(t, x, w*h*sizeof*t);
	float sa = -180*a/M_PI;
	//if (a > 0) sa = 180+sa;
	cast_horizontal_shadows(t, w, h, sa);

	// TODO: find reflectors and simulate them
	// ...

	// project stuff
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		if (!isfinite(t[j*w+i])) continue;
		float p = w/2 + cos(a) * (i-w/2) - sin(a) * x[j*w+i];
		int ip = floor(p);
		float q = hypot(1, hypot(x[j*w+i+1]-x[j*w+i], x[j*w+i+w]-x[j*w+i]));
		q = 1;
		float fp = (p - ip);
		//fprintf(stderr, "(!) p=%g ip=%d fp=%g\n", p, ip, fp);
		assert(fp >= 0);
		assert(fp < 1);
		//float r = 1; // may be random number
		//float r = random_uniform();
		//float r0 = random_uniform();
		if (ip >= 0 && ip < w-1)
		{
			y[j*w+ip+0] += q * (1-fp) * pow(random_uniform(),2);
			y[j*w+ip+1] += q * fp     * pow(random_uniform(),2);
		}
		else
			if (j == 10)
				fprintf(stderr, "bad ip=%d\n", ip);
	}
	free(t);
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
