typedef float (*getpixel_operator)(float*,int,int,int,int);
typedef void (*setpixel_operator)(float*,int,int,int,int,float);

static float getpixel_127(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 127;
	return x[i + j*w];
}

static void setpixel_insideP(float *x, int w, int h, int i, int j, float v)
{
	if (i >= 0 && j >= 0 && i < w && j < h)
		x[w*j+i] = v;
}

static void setpixel_tsum_insideP(float *x, int w, int h, int i, int j, float v)
{
	if (i >= 0 && j >= 0 && i < w && j < h)
	{
		float g = x[w*j+i] + v;
		x[w*j+i] = g > 255 ? 255 : g;
	}
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

// Floyd-Steinberg
void dither(float *x, int w, int h)
{
	getpixel_operator get = getpixel_127;
	setpixel_operator set = setpixel_insideP;
	setpixel_operator add = setpixel_tsum_insideP;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float old = get(x, w, h, i, j);
		float new = dither_value(old);
		float err = old - new;
		set(x, w, h, i, j, new);
		add(x, w, h, i+1, j+0, 7*err/16);
		add(x, w, h, i-1, j+1, 3*err/16);
		add(x, w, h, i+0, j+1, 5*err/16);
		add(x, w, h, i+1, j+1, 1*err/16);
	}
}

void dither_sep(float *x, int w, int h, int pd)
{
	for (int i = 0; i < pd; i++)
		dither(x + i*w*h, w, h);
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [gray [binary]]\n", *v);
		//                          0  1     2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	dither_sep(x, w, h, pd);


	iio_save_image_float_split(filename_out, x, w, h, pd);

	free(x);
	return 0;
}
