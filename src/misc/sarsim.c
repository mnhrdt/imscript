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

static float random_speckle(void)
{
	//return 1;
	//return random_uniform();
	return pow(random_uniform(),2);
	//return pow(random_normal(),2);
	//return fabs(random_normal());
	//return random_uniform()*random_uniform();
}

static void radar_sim_horizontal(
		float *y,      // output: radar simulation
		float *x,      // input: height raster in pixel units
		int w,         // raster width
		int h,         // raster height
		float a        // angle pi=left pi/2=top 0=right
		)
{
	// TODO: rewrite this whole function as a random sampling across the
	// whole image, computing the surface element at each point and acting
	// accordingly (after taking into account the shadows, of course)

	fprintf(stderr, "radarsim a=%g\n", a);

	// initialize output canvas to zero
	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	// remove invisible parts (preprocess with shadowcast)
	float *t = malloc(w*h*sizeof*t); // the shadow mask
	memcpy(t, x, w*h*sizeof*t);
	float sa = -180*a/M_PI;
	cast_horizontal_shadows(t, w, h, sa);
	fprintf(stderr, "a=%g sa=%g\n", a, sa);

	// TODO: find reflectors and simulate them
	//if (a > 0) { // left to right
	//	for (int j = 0; j < h-1; j++)
	//	for (int i = 0; i < w-1; i++)
	//	{
	//		if (!isfinite(t[j*w+i])) continue;
	//		if (!isfinite(t[j*w+i+1])) continue;
	//		float x0 = x[j*w+i+0];
	//		float x1 = x[j*w+i+1];
	//		if (x1 - x0 > 2) // tall enough vertical wall
	//		{
	//			float p = w/2 + cos(a) * (i-w/2) - sin(a) * x0;
	//			int ip = floor(p);
	//			if (ip >= 0 && ip < w-1)
	//				y[j*w+ip+0] += pow(x1-x0, 1);
	//		}
	//	}
	//}
	//if (a < 0) { // right to left
	//}


	// TODO: find vertical walls and sample their points
	if (a > 0) { // left to right
		for (int j = 0; j < h-1; j++)
		for (int i = 0; i < w-1; i++)
		{
			if (!isfinite(t[j*w+i])) continue;
			if (!isfinite(t[j*w+i+1])) continue;
			float x0 = x[j*w+i+0];
			float x1 = x[j*w+i+1];
			if (x1 - x0 > 2) // tall enough vertical wall
			for (int k = 0; k < lrint(x1-x0); k++)
			{
				float p = w/2 + cos(a)*(i-w/2) - sin(a)*(x0+k);
				int ip = floor(p);
				if (ip >= 0 && ip < w)
					y[j*w+ip+0] += 1*random_speckle();
			}
		}
	}
	if (a < 0) { // right to left
		for (int j = 0; j < h-1; j++)
		for (int i = w-1; i > 1; i--)
		{
			if (!isfinite(t[j*w+i])) continue;
			if (!isfinite(t[j*w+i-1])) continue;
			float x0 = x[j*w+i+0];
			float x1 = x[j*w+i-1];
			if (x1 - x0 > 2) // tall enough vertical wall
			for (int k = 0; k < lrint(x1-x0); k++)
			{
				float p = w/2 + cos(a)*(i-w/2) - sin(a)*(x0+k);
				int ip = floor(p);
				if (ip >= 0 && ip < w)
					y[j*w+ip+0] += 1*random_speckle();
			}
		}
	}

	// project stuff
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		if (!isfinite(t[j*w+i])) continue;
		float p = w/2 + cos(a) * (i-w/2) - sin(a) * x[j*w+i];
		int ip = floor(p);
		float q = hypot(1, hypot(x[j*w+i+1]-x[j*w+i], x[j*w+i+w]-x[j*w+i]));
		q = fmin(q, 3);
		float fp = (p - ip);
		//fprintf(stderr, "(!) p=%g ip=%d fp=%g\n", p, ip, fp);
		assert(fp >= 0);
		assert(fp < 1);
		//float r = 1; // may be random number
		//float r = random_uniform();
		//float r0 = random_uniform();
		if (ip >= 0 && ip < w-1)
		{
			y[j*w+ip+0] += q * (1-fp) * random_speckle();
			y[j*w+ip+1] += q * fp     * random_speckle();
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
