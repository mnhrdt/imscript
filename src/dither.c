

typedef float (*getpixel_operator)(float*,int,int,int,int);
typedef void (*setpixel_operator)(float*,int,int,int,int,float);

static float getpixel_127(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

static void setpixel_insideP(float *x, int w, int h, int i, int j, float v)
{
	if (i >= 0 && j >= 0 && i < w && j < h)
		x[w*j+i] = v;
}

static void setpixel_sum_insideP(float *x, int w, int h, int i, int j, float v)
{
	if (i >= 0 && j >= 0 && i < w && j < h)
		x[w*j+i] += v;
}

static float dither_value(float x)
{
	return x < 127.5 ? 0 : 255;
}

void dither(float *x, int w, int h)
{
	getpixel_operator get = getpixel_127;
	setpixel_operator set = setpixel_insideP;
	setpixel_operator add = setpixel_sum_insideP;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float o = get(x, w, h, i, j);
		float n = dither_value(o);
		set(x, w, h, i, j, n);
		float e = o - n;
		add(x, w, h, i+1, j+0, 7*e/16);
		add(x, w, h, i-1, j+1, 3*e/16);
		add(x, w, h, i+0, j+1, 5*e/16);
		add(x, w, h, i+1, j+1, 1*e/16);
	}
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c > 1) {
		fprintf(stderr, "usage:\n\t%s < gray > binary\n", *v);
		return 1;
	}
	char *filename_in = "-";
	char *filename_out = "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	dither(x, w, h);

	iio_save_image_float(filename_out, x, w, h);

	free(x);
	return 0;
}
