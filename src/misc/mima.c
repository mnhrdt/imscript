// set to zero, except positive local maxima and negative local minima

#include <math.h>
#include <stdbool.h>

float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return NAN;
	return x[j*w+i];
}

void mima(float *y, float *x, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		bool rmax = true, rmin = true;
		for (int dj = -1; dj <= 1; dj++)
		for (int di = -1; di <= 1; di++)
		{
			if (!di && !dj) continue;
			float xij = getpixel_nan(x, w, h, i   , j   );
			float xIJ = getpixel_nan(x, w, h, i+di, j+dj);
			if (xij <= 0 || xij <= xIJ) rmax = false;
			if (xij >= 0 || xij >= xIJ) rmin = false;
		}
		y[j*w+i] = (rmax || rmin) ? x[j*w+i] : 0;
	}
}

void mima_separable(float *y, float *x, int w, int h, int pd)
{
	for (int l = 0; l < pd; l++)
		mima(y + w*h*l, x + w*h*l, w, h);
}

#include <stdlib.h>
#include <stdio.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3)
		return fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
	//                                         0  1   2
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = malloc(w*h*pd*sizeof*y);

	mima_separable(y, x, w, h, pd);

	iio_write_image_float_split(filename_out, y, w, h, pd);

	free(x);
	free(y);
	return 0;
}
