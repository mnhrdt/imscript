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

// boustrophedonic Floyd-Steinberg
void dither_b(float *x, int w, int h)
{
	getpixel_operator get = getpixel_127;
	setpixel_operator set = setpixel_insideP;
	setpixel_operator add = setpixel_tsum_insideP;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = (j%2) ? w - i - 1 : i;
		float old = get(x, w, h, ii, j);
		float new = dither_value(old);
		float err = old - new;
		set(x, w, h, ii, j, new);
		if (j%2) {
			add(x, w, h, ii-1, j+0, 7*err/16);
			add(x, w, h, ii+1, j+1, 3*err/16);
			add(x, w, h, ii-0, j+1, 5*err/16);
			add(x, w, h, ii-1, j+1, 1*err/16);
		} else {
			add(x, w, h, ii+1, j+0, 7*err/16);
			add(x, w, h, ii-1, j+1, 3*err/16);
			add(x, w, h, ii+0, j+1, 5*err/16);
			add(x, w, h, ii+1, j+1, 1*err/16);
		}
	}
}

// Floyd-Steinberg "adaptive", by github user zvezdochiot
void dither_a(float *x, int w, int h)
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
		float val = get(x, w, h, i-1, j+0);
		val = (val > new) ? (val - new) : (new - val);
		float errk1 = 6*val + 1;
		val = get(x, w, h, i+1, j+1);
		val = (val > new) ? (val - new) : (new - val);
		float errk2 = 2*val + 1;
		val = get(x, w, h, i+0, j+1);
		val = (val > new) ? (val - new) : (new - val);
		float errk3 = 4*val + 1;
		float errk4 = 1;
		float errks = errk1 + errk2 + errk3 + errk4;
		add(x, w, h, i-1, j+0, err*errk1/errks);
		add(x, w, h, i+1, j+1, err*errk2/errks);
		add(x, w, h, i-0, j+1, err*errk3/errks);
		add(x, w, h, i-1, j+1, err*errk4/errks);
		set(x, w, h, i, j, new);
	}
}


void dither_sep(float *x, int w, int h, int pd, _Bool b)
{
	for (int i = 0; i < pd; i++)
		if (b)
			dither_b(x + i*w*h, w, h);
		else
			dither(x + i*w*h, w, h);
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"

static char *help_string_name     = "dither";
static char *help_string_version  = "dither 1.0\n\nWritten by eml";
static char *help_string_oneliner = "binarize a grayscale image by dithering";
static char *help_string_usage    = "usage:\n\t"
"dither [-b] [in [out]]";
static char *help_string_long     =
"Filter that binarizes a grayscale image by the Floyd-Steinberg algorithm.\n"
"\n"
"The input image is assumed to take values between 0=black and 255=white.\n"
"The output image is binary on the same range.  For color images,\n"
"the algorithm is applied independently to each color channel.\n"
"\n"
"Usage: dither in.png out.png\n"
"   or: dither in.png > out.npy\n"
"   or: cat in.png | dither > out.npy\n"
"\n"
"Options:\n"
" -b\t\tboustrophedonic traversal order\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" dither lena.png lena.pbm             Floyd-Steinberg dithering\n"
" blur C 1 lena.png | qauto | dither   contrast-enhanced dithering\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
int main_dither(int c, char *v[])
{
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	_Bool option_b = pick_option(&c, &v, "b", 0);
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [gray [binary]]\n", *v);
		//                          0  1     2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	dither_sep(x, w, h, pd, option_b);


	iio_write_image_float_split(filename_out, x, w, h, pd);

	free(x);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_dither(c, v); }
#endif
