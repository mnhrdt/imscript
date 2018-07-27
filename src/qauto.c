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

static void get_rminmax(float *rmin, float *rmax, float *x, int n, float rb)
{
	float *tx = xmalloc(n*sizeof*tx);
	int N = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(x[i]))
			tx[N++] = x[i];
	int irb = round(rb);
	if (N < 1) {
		fprintf(stderr, "too many NANs (rb N) = %g %d", rb, N);
		abort();
	}
	qsort(tx, N, sizeof*tx, compare_floats);
	*rmin = tx[irb];
	*rmax = tx[N-1-irb];
	free(tx);
}

static void get_avgstd(float *out_avg, float *out_std, float *x, int n)
{
	long double avg = 0;
	for (int i = 0; i < n; i++)
		avg += x[i];
	avg /= n;
	*out_avg = avg;

	long double var = 0;
	for (int i = 0; i < n; i++)
		var += (x[i] - avg) * (x[i] - avg);
	var /= n;
	*out_std = sqrt(var);
}

// adjust contrast by allowing a fixed percentile of top and bottom saturation
static void simplest_color_balance(float *y, float *x, int w, int h, float p)
{
	float rmin, rmax;
	get_rminmax(&rmin, &rmax, x, w*h, w*h*(p/100));
	if (global_verbose_flag)
		fprintf(stderr, "qauto: rminmax = %g %g\n", rmin, rmax);

	for (int i = 0; i < w*h; i++)
		y[i] = 255 * (x[i] - rmin) / (rmax - rmin);
}

// adjust contrast by normalizing avg=127 std=s
static void adjust_avgstd(float *y, float *x, int w, int h, float s)
{
	float avg, std;
	get_avgstd(&avg, &std, x, w*h);
	fprintf(stderr, "qauto: avgstd = %g %g\n", avg, std);

	for (int i = 0; i < w*h; i++)
		y[i] = 127 + s * (x[i] - avg) / std;
}

static void qauto_grey(float *y, float *x, int w, int h, float parameter)
{
	if (parameter >= 0)
		simplest_color_balance(y, x, w, h, parameter);
	else
		adjust_avgstd(y, x, w, h, -parameter);
}

static void quantize_to_byte_values_inplace(float *x, int n)
{
	for (int i = 0; i < n; i++)
	{
		float g = round(x[i]);
		int ig = g;
		if (ig < 0) ig = 0;
		if (ig > 255) ig = 255;
		x[i] = ig;
	}
}

static void qauto(float *y, float *x, int w, int h, int pd,
	bool independent_channels, float parameter, bool do_not_quantize)
{
	if (independent_channels)
		for (int l = 0; l < pd; l++)
			qauto_grey(y + w*h*l, x + w*h*l, w, h, parameter);
	else
		qauto_grey(y, x, w, h*pd, parameter);

	if (!do_not_quantize)
		quantize_to_byte_values_inplace(y, w*h*pd);
}


static char *help_string_name     = "qauto";
static char *help_string_version  = "qauto 1.0\n\nWritten by eml";
static char *help_string_oneliner = "quantize an image into [0,255]";
static char *help_string_usage    = "usage:\n\t"
"qauto [-p 5] [-i] [-f]  [in [out]]";
static char *help_string_long     =
"Qauto quantizes an image into [0,255].\n"
"\n"
"The image is trasformed by an affine contrast change I -> a*I + b\n"
"and then the colors are saturated and quantized into [0,255].\n"
"The parameters (a,b) of the contrast change are computed to statisfy\n"
"certain conditions. By default, they are chosen so that 5% of the pixels\n"
"are saturated.\n"
"\n"
"Usage: qauto in.tiff out.png\n"
"   or: qauto in.tiff > out.pnm\n"
"   or: cat in.tiff | qauto > out.pnm\n"
"\n"
"Options:\n"
" -p X\t\tuse a percentile of X% (default X=5)\n"
" -p -X\t\tset the mean to 127 and the standard deviation to X\n"
" -f\t\tdo not quantize the output, only rescale the values\n"
" -i\t\ttreat each pixel dimension independently\n"
" -v\t\tverbose mode (print details of the transformation)\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" qauto in.tiff out.png          Quantize an image by simplest color balance.\n"
" qauto -i in.png out.png        Remove color biases\n"
"\n"
"Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"    // function to extract hyphenated command line options
int main_qauto(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// extract named options
	float parameter = atof(pick_option(&c, &v, "p", "5"));
	bool independent = pick_option(&c, &v, "i", 0);
	bool dontquantize = pick_option(&c, &v, "f", 0);
	global_verbose_flag = !!pick_option(&c, &v, "v", 0);

	// get positional arguments
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t"
				"%s [in [out]] {-p 5} {-i} {-f}\n", *v);
		//                0  1   2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(in, &w, &h, &pd);

	// allocate space for output image
	float *y = xmalloc(w*h*pd*sizeof*y);

	// run the algorithm
	qauto(y, x, w, h, pd, independent, parameter, dontquantize);

	// write result and exit
	iio_write_image_float_split(out, y, w, h, pd);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_qauto(c, v); }
#endif
