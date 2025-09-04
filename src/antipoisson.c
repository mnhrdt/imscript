// poisson pixel sampler with a given exclusion area

#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "random.c"


static int good_modulus(int n, int p)
{
	int r = n % p;
	r = r < 0 ? r + p : r;
	assert(r >= 0);
	assert(r < p);
	return r;
}

static void antipoisson(float *y, float *x, int w, int h, int n, int t)
{
	// start with a black canvas
	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	// build the list of masked pixels
	int m = 0;  // number of pixels in the mask
	for (int i = 0; i < w*h; i++)
		if (x[i])
			m += 1;
	int M[m][2];  // list of masked pixel positions
	int m2 = 0;
	for (int i = 0; i < w*h; i++)
		if (x[i])
		{
			M[m2][0] = i % w;
			M[m2][1] = i / w;
			m2 += 1;
		}

	// for each of the n requested points, perform t trial samples
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < t; j++)
		{
			int p = randombounds(0, w-1);
			int q = randombounds(0, h-1);
			if (!y[q*w+p])
			{
				for (int k = 0; k < m; k++)
				{
					int pp = good_modulus(p + M[k][0], w);
					int qq = good_modulus(q + M[k][1], h);
					y[qq*w+pp] += 1;
				}
				y[q*w+p] = 255;
				break; // found a point, move to next
			}
		}
	}
}



static char *help_string_name     = "antipoisson";
static char *help_string_version  = "antipoisson 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "Poisson sampler with exclusion region";
static char *help_string_usage    = "usage:\n\t"
"antipoisson [in_mask [out_pix]]";
static char *help_string_long     =
"Antipoisson give a binary image of dots sampled by avoiding a local mask "
"\n"
"Usage: antipoisson in.npy out.png\n"
"   or: antipoisson in.npy > out.pnm\n"
"   or: cat in.npy | antipoisson > out.pnm\n"
"\n"
"Options:\n"
" -n N\t\tnumber of points (default 100)\n"
" -t T\t\tnumber of failed trials before giving out (default 100)\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" antipoisson mask.png dots.png.\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"    // function to extract hyphenated command line options
int main_antipoisson(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// extract named options
	int n = atoi(pick_option(&c, &v, "n", "100"));
	int t = atoi(pick_option(&c, &v, "t", "100"));

	// get positional arguments
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t"
				"%s [in [out]]\n", *v);
		//                0  1   2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	// read input mask
	int w, h;
	float *x = iio_read_image_float(in, &w, &h);

	// allocate space for output image
	float *y = malloc(w*h*sizeof*y);

	// run the algorithm
	antipoisson(y, x, w, h, n, t);

	// write result and exit
	iio_write_image_float_vec(out, y, w, h, 1);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_antipoisson(c, v); }
#endif
