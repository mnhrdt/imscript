// monodisp: remove non-monotonic points in a disparity map

#include <math.h> // NAN

void monodisp(float *y, float *x, int w, int h, char *criterion)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float *d = x + w*j;
		if (i > 0 && d[i] - d[i-1] <= -1)
		{
			if (*criterion == 'l' || *criterion == 'b')
				y[j*w+i-1] = NAN;
			if (*criterion == 'r' || *criterion == 'b')
				y[j*w+i] = NAN;
		}
		else
			y[j*w+i] = d[i];
	}
}

#define MAIN_MONODISP

#ifdef MAIN_MONODISP
#include "iio.h"
#include <stdio.h>
#include <stdlib.h>
int main(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 3 && c != 4)
		return fprintf(stderr,
			"usage:\n\t%s {left|right|both} [in [out]]\n", *v);
			//          0  1                 2   3
	char *criterion    = v[1];
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// load input image, allocate space for output image
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *y = malloc(w*h*sizeof*y);

	// run the algorithm
	monodisp(y, x, w, h, criterion);

	// save and quit
	iio_write_image_float(filename_out, y, w, h);
	return 0;
}
#endif//MAIN_MONODISP
