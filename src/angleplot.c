#include <stdio.h>
#include <stdlib.h>

static void angleplot(float *yy, int width, float scale,
		float *xx, float w, float h)
{
	float (*y)[width] = (void*)yy;
	float (*x)[2] = (void*)xx;

	// background
	for (int i = 0; i < width*width; i++)
		yy[i] = 0;

	// axes
	for (int i = 0; i < width; i++)
		y[width/2][i] = y[i][width/2] = 1;

	// data points
	float alpha = width/(2*scale);
	float beta = width/2.0;
	for (int i = 0; i < w * h; i++)
	{
		int q[2];
		for (int j = 0; j < 2; j++) {
			q[j] = alpha * x[i][j] + beta;
			if (q[j] < 0) q[j] = 0;
			if (q[j] >= width) q[j] = width-1;
		}
		y[q[1]][q[0]] += 1;
	}
}

#ifndef OMIT_MAIN



#include "iio.h"

#include "fail.c"
#include "xmalloc.c"

int main(int c, char *v[])
{
	if (c < 3 || c > 5) {
		fprintf(stderr, "usage:\n\t%s width scale [flow [view]]\n", *v);
		//                          0 1     2      3     4
		return EXIT_FAILURE;
	}
	int width = atoi(v[1]);
	float scale = atof(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";
	int w, h, pd;
	float *f = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (pd != 2) fail("need 2D-valued input");
	float *o = xmalloc(width * width * sizeof*o);
	angleplot(o, width, scale, f, w, h);
	iio_save_image_float_vec(filename_out, o, width, width, 1);
	return EXIT_SUCCESS;
}

#endif//OMIT_MAIN
