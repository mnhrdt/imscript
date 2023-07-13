#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

static bool numericP(float *x, int d)
{
	for (int i = 0; i < d; i++)
		if (isnan(x[i]))
			return false;
	return true;
}

static bool nonzeroP(float *x, int d)
{
	for (int i = 0; i < d; i++)
		if (x[i] != 0)
			return true;
	return false;
}

void columnize(FILE *f, float *x, int w, int h, int pd, bool pz, bool pv)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	if (numericP(x + (j*w + i)*pd, pd))
	if (nonzeroP(x + (j*w + i)*pd, pd) || pz)
	{

		fprintf(f, "%d\t%d", i, j);
		if (pv)
		for (int k = 0; k < pd; k++)
			fprintf(f, "\t%g", x[(j*w+i)*pd + k]);
		fprintf(f, "\n");
	}
}


static char *help_string_name     = "columnize";
static char *help_string_version  = "columnize 1.0\n\nWritten by eml";
static char *help_string_oneliner = "describe an image as a multi-column file";
static char *help_string_usage    = "usage:\n\t"
"columnize {-z} [in.img [out.txt]]";
static char *help_string_long     =
"Columnize prints the positions and values of all pixels of an image\n"
"\n"
"By default, zero-valued pixels are not printed.  If you want to print\n"
"them, use the -z option.  Nan pixels are never printed.\n"
"\n"
"Usage: columnize in.png out.txt\n"
"   or: columnize in.png > out.txt\n"
"   or: cat in.png | columnize > out.txt\n"
"\n"
"Options:\n"
" -z\t\tprint also zero-valued pixels\n"
" -n\t\tdo not print the pixel values, only the positions\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" columnize in.png | wc -l       Count the non-zero pixels of an image.\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"    // function to extract hyphenated command line options
#include "xfopen.c"     // open a file wrapper
int main_columnize(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// extract named options
	bool printzero = pick_option(&c, &v, "z", 0);
	bool printvals = !pick_option(&c, &v, "n", 0);

	// get positional arguments
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t"
				"%s [in [out]] {-z} {-n}\n", *v);
		//                0  1   2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);

	// run the algorithm
	FILE *f = xfopen(out, "w");
	columnize(f, x, w, h, pd, printzero, printvals);

	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_columnize(c, v); }
#endif
