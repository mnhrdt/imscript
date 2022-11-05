#include <stdio.h> // only for "fprintf"
#include <stdlib.h> // only for "free"
#include <stdint.h>
#include "iio.h"

// read an image from a file or from stdin, and write the named file
int main_iion(int c, char *v[])
{
	if (c != 3 && c != 2)
		return fprintf(stderr, "usage:\n\t%s [infile] outfile\n", *v);
	//                                         0 1        2
	char *filename_in = c < 3 ? "-" : v[1];
	char *filename_out = v[c - 1];
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (!x)
		return fprintf(stderr, "cannot image \"%s\"\n", filename_in);
	iio_write_image_float_vec(filename_out, x, w, h, pd);
	free(x);
	return 0;
}
#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_iion(c, v); }
#endif
