#include <stdlib.h> // atof
#include <stdio.h>  // fprintf
#include <math.h>   // floor
#include "iio.h"    // iio_*

static char *help_string_name     = "qeasy";
static char *help_string_version  = "qeasy 1.0\n\nWritten by eml";
static char *help_string_oneliner = "quantize an image into [0,255]";
static char *help_string_usage    = "usage:\n\t"
"qeasy black white [in [out]]";
static char *help_string_long     =
"Qeasy quantizes an image into [0,255] using the selected range.\n"
"\n"
"The image is trasformed by an affine contrast change I -> a*I + b\n"
"and then the colors are saturated and quantized into [0,255].\n"
"The parameters (a,b) of the contrast change are computed to statisfy\n"
"certain conditions. By default, they are chosen so that 5% of the pixels\n"
"are saturated.\n"
"\n"
"Usage: qeasy black white in.tiff out.png\n"
"   or: qeasy black white in.tiff > out.pnm\n"
"   or: cat in.tiff | qeasy black white > out.pnm\n"
"\n"
"Options:\n"
" -f\t\tdo not quantize the output, only rescale the values\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" qauto in.tiff out.png          Quantize an image by simplest color balance.\n"
" qauto -i in.png out.png        Remove color biases\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
int main_qeasy(int c, char *v[])
{
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s black white  [in [out]]\n", *v);
		//                         0 1     2       3   4
		return 1;
	}
	float black = atof(v[1]);
	float white = atof(v[2]);
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	uint8_t *y = (void*)x;
	for (int i = 0; i < w*h*pd; i++) {
		float g = x[i];
		g = floor(255 * (g - black)/(white - black));
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		y[i] = g;
	}
	iio_write_image_uint8_vec(filename_out, y, w, h, pd);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_qeasy(c, v); }
#endif
