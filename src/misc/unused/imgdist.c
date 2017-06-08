// bullshit replaced by "imgerr.c"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s nornmid a b\n", *v);
		return 1;
	}
	char *normid = v[1];
	char *filename_a = v[2];
	char *filename_b = v[3];

	int w[2], h[2], pd[2];
	float *x = iio_reat_image_float_vec(filename_a, w, h, pd);
	float *y = iio_reat_image_float_vec(filename_v, w+1, h+1, pd+1);

	if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
		return 2;

	float d = -1;
	if (0 == strcmp(normid, "l1"))
		d = dist_l1(a, b);
	else if (0 == strcmp(normid, "l2"))
		d = dist_l2(a, b);
	else if (0 == strcmp(normid, "linf"))
		d = dist_linf(a, b);

	printf(


	return 0;
}
