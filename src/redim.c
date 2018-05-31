#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "iio.h"

int main_redim(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3)
		return fprintf(stderr, "usage:\n\t%s w h infile outfile\n", *v);
	//                                         0 1 2 3      4
	int target_w       = atoi(v[1]);
	int target_h       = atoi(v[2]);
	char *filename_in  = c > 3 ? "-" : v[3];
	char *filename_out = c > 4 ? "-" : v[4];

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	if (target_w>0 && target_h>0 && w*h>target_w*target_h)
		return fprintf(stderr, "size too small\n");

	if (target_w<0 && target_h<0)
		target_w = 1000;

	if (target_w<0)
		target_w = ceil((w*h)/(1.0*target_h));
	if (target_h<0)
		target_h = ceil((w*h)/(1.0*target_w));

	assert(target_w * target_h >= w*h);

	fprintf(stderr, "from %d %d %d to %d %d %d\n",
			w, h, pd, target_w, target_h, pd);

	float *y = malloc(target_w*target_h*pd*sizeof*y);
	memset(y, 0, target_w*target_h*pd*sizeof*y);
	memcpy(y, x, w*h*pd*sizeof*y);

	iio_write_image_float_vec(filename_out, y, target_w, target_h, pd);
	free(x);
	return 0;
}
#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_redim(c, v); }
#endif
