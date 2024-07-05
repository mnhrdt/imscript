// read an image from a file or from stdin, and write the named file

#include <stdio.h>   // fprintf
#include <stdlib.h>  // free
#include <stdint.h>  // uint8_t, uint16_t
#include <string.h>  // strcmp
#include "iio.h"     // iio_read_image_T_vec, iio_write_image_T_vec

typedef uint8_t uint8;
typedef uint16_t uint16;

#define X(t) int main_iion_ ## t(int c, char *v[])                          \
{                                                                           \
	if (c != 3 && c != 2)                                               \
		return fprintf(stderr, "usage:\n\t%s [in] out\n", *v);      \
	char *fname_in = c < 3 ? "-" : v[1];                                \
	char *fname_out = v[c - 1];                                         \
	int w, h, pd;                                                       \
	t *x = iio_read_image_ ## t ## _vec(fname_in, &w, &h, &pd);         \
	if (!x)                                                             \
		return fprintf(stderr, "cannot image \"%s\"\n", fname_in);  \
	iio_write_image_ ## t ## _vec(fname_out, x, w, h, pd);              \
	free(x);                                                            \
	return 0;                                                           \
}
X(float)
X(double)
X(uint8)
X(uint16)
#undef X


static char *help_string_name     = "iion";
static char *help_string_version  = "iion 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "read and write a named image using iio";
static char *help_string_usage    = "usage:\n\t"
"iion [-t {float,double,uint8,uint16}] [in] out";
static char *help_string_long     =
"Iion reads an image from a file or from stdin, and writes the named file.\n"
"\n"
"This is a \"no-op\" filter for images.  The sole application is to convert\n"
"between image file formats according to the semantics of iio(3).\n"
"The default type is float, but this can be changed by options.\n"
"\n"
"Usage: iion in out.png\n"
"   or: cat in | iion out.png\n"
"\n"
"Options:\n"
" -t T\t\tuse T as intermediate type (float,double,uint8,uint16)\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" iion http://path/to/img x.png       Download an image as png\n"
" iion -t uint16 x.jpg x.tiff         Convert jpeg to 16 bit tiff.\n"
" cat img.npy | iion img.png          Quantize float data to 8 bits\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"    // function to extract hyphenated command line options
int main_iion(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);

	char *t = pick_option(&c, &v, "t", "float");

#define X(x) if (0 == strcmp(t, #x)) return main_iion_ ## x(c, v)
	X(float);
	X(double);
	X(uint8);
	X(uint16);
#undef X

	return fprintf(stderr, "%s: unrecognized data type \"%s\"\n", *v, t);
}


#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_iion(c, v); }
#endif
