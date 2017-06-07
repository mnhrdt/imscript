

typedef float (*getpixel_operator)(float*,int,int,int,int);
typedef void (*setpixel_operator)(float*,int,int,int,int,float);

static float getpixel_0(float *x, int w, int h, int i, int j)
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

static float dither_value(float x, float lev, float top)
{
	return x < lev ? 0 : top;
}

// Floyd-Steinberg
void dither(float *x, int w, int h, float lev, float top)
{
	getpixel_operator get = getpixel_0;
	setpixel_operator set = setpixel_insideP;
	setpixel_operator add = setpixel_sum_insideP;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float old = get(x, w, h, i, j);
		float new = dither_value(old, lev, top);
		float err = old - new;
		set(x, w, h, i, j, new);
		add(x, w, h, i+1, j+0, 7*err/16);
		add(x, w, h, i-1, j+1, 3*err/16);
		add(x, w, h, i+0, j+1, 5*err/16);
		add(x, w, h, i+1, j+1, 1*err/16);
	}
}

void dither_sep(float *x, int w, int h, int pd, float lev, float top)
{
	for (int i = 0; i < pd; i++)
		dither(x + i*w*h, w, h, lev, top);
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *lev_string = pick_option(&c, &v, "l", "127.5");
	char *top_string = pick_option(&c, &v, "t", "255");
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t"
				"%s [-l 127] [-t 255] [gray [binary]]\n", *v);
		//                0                    1     2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";
	float lev = atof(lev_string);
	float top = atof(top_string);

	fprintf(stderr, "lev = %g\n", lev);
	fprintf(stderr, "top = %g\n", top);

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	dither_sep(x, w, h, pd, lev, top);

	iio_write_image_float_split(filename_out, x, w, h, pd);

	free(x);
	return 0;
}
