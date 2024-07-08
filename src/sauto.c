#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define xmalloc malloc
static int global_verbose_flag = 0;

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static float clip(float x, float a, float b)
{
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static void sauto(uint8_t *y, float *x, int w, int h, float p)
{
	int n = 0; // number of non-nan samples
	float *t = xmalloc(w*h*sizeof*t); // table of numeric samples (to sort)
	for (int i = 0; i < w*h; i++)
		if (!isnan(x[i]))
			t[n++] = fabs(x[i]);
	qsort(t, n, sizeof*t, compare_floats);
	float s = 0; // saturation quantile
	if (p >= 0) // p is a percentage
	{
		//assert(p <= 50);
		int i = n - 1 - p*n/100;
		if (i < 0) i = 0;
		if (i >= n) i = n;
		//assert(i >= 0);
		//assert(i < n);
		s = t[i];
	} else { // -p is a number of pixels
		int i = n + p - 1;
		if (i < 0) i = 0;
		if (i >= n) i = n - 1;
		fprintf(stderr, "n=%d p=%g i=%d\n", n, p, i);
		s = t[i];
	}
	if (global_verbose_flag) fprintf(stderr, "sauto: s = %g\n", s);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		uint8_t rgb[3] = {100, 100, 100}; // color for NAN
		float v = x[j*w+i];
		if (!isnan(v))
		{
			float r = 1 - clip(v/s, 0, 1);
			float g = 1 - clip(fabs(v/s), 0, 1);
			float b = 1 + clip(v/s, -1, 0);
			rgb[0] = 255*r;
			rgb[1] = 255*g;
			rgb[2] = 255*b;
		}
		for (int k = 0; k < 3; k++)
			*y++ = rgb[k];
	}
}


static char *help_string_name     = "sauto";
static char *help_string_version  = "sauto 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "render a signed image with color palette";
static char *help_string_usage    = "usage:\n\t"
"sauto [-p 1] [in [out]]";
static char *help_string_long     =
"Sauto applies a blue-white-red palette to a signed image.\n"
"\n"
"The image is trasformed by an affine contrast change I -> a*I + b\n"
"and then the colors are saturated and quantized into [0,255].\n"
"The parameters (a,b) of the contrast change are computed to statisfy\n"
"certain conditions. By default, they are chosen so that 5% of the pixels\n"
"are saturated.\n"
"\n"
"Usage: sauto in.npy out.png\n"
"   or: qauto in.npy > out.pnm\n"
"   or: cat in.npy | qauto > out.pnm\n"
"\n"
"Options:\n"
" -p X\t\tsaturate a percentile of X% (default X=1)\n"
" -v\t\tverbose mode (print details of the transformation)\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" qauto in.npy out.png          Quantize an image by simplest color balance.\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"    // function to extract hyphenated command line options
int main_sauto(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// extract named options
	float p = atof(pick_option(&c, &v, "p", "1"));
	global_verbose_flag = !!pick_option(&c, &v, "v", 0);

	// get positional arguments
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t"
				"%s [in [out]] {-p 1}\n", *v);
		//                0  1   2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	// read input image
	int w, h;
	float *x = iio_read_image_float(in, &w, &h);

	// allocate space for output image
	uint8_t *y = xmalloc(3*w*h);

	// run the algorithm
	sauto(y, x, w, h, p);

	// write result and exit
	iio_write_image_uint8_vec(out, y, w, h, 3);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_sauto(c, v); }
#endif
