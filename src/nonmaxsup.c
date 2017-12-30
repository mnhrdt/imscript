// non-maximum suppression

#include <math.h>

#include "smapa.h"
SMART_PARAMETER(GRADFAC,1)

void non_maximum_suppression(
		float *o,      // output (copy of x with non-maxima set to 0)
		float *x,      // input scalar field
		float *g,      // input vector field (to look for neighbors)
		int w,         // image width
		int h          // image height
		)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float v0 = x[j*w+i];
		float v1 = -INFINITY;
		float v2 = -INFINITY;
		float gx = g[2*(j*w+i)+0];
		float gy = g[2*(j*w+i)+1];
		float ngx = gx / hypot(gx, gy);
		float ngy = gy / hypot(gx, gy);
		float f = GRADFAC();
		{
			int ii = lrint(i + f * ngx);
			int jj = lrint(j + f * ngy);
			if (ii >= 0 && jj >= 0 && ii < w && jj < h)
				v1 = x[jj*w+ii];
		}
		{
			int ii = lrint(i - f * ngx);
			int jj = lrint(j - f * ngy);
			if (ii >= 0 && jj >= 0 && ii < w && jj < h)
				v2 = x[jj*w+ii];
		}
		o[j*w+i] = (v0 > v1 && v0 > v2) ? v0 : 0;
	}
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 4)
		return fprintf(stderr, "usage:\n\t%s scal vect outscal\n", *v);
		//                                 0 1    2    3
	char *filename_scal = v[1];
	char *filename_vect = v[2];
	char *filename_out  = v[3];

	// read input images
	int w[2], h[2], pd[2];
	float *x = iio_read_image_float_vec(filename_scal, w+0, h+0, pd+0);
	float *y = iio_read_image_float_vec(filename_vect, w+1, h+1, pd+1);
	if (w[0] != w[1] || h[0] != h[1] || pd[0] != 1 || pd[1] != 2)
		return fprintf(stderr, "bad dims (%d %d %d),(%d %d %d)\n",
				w[0], h[0], pd[0], w[1], h[1], pd[1]);

	// allocate space for output
	float *z = malloc(*w * *h * sizeof*z);

	// perform computation
	non_maximum_suppression(z, x, y, *w, *h);

	// write output image
	iio_write_image_float_vec(filename_out, z, *w, *h, 1);

	// cleanup and exit
	free(x);
	free(y);
	free(z);
	return 0;
}
