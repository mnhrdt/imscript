static void compute_intimg_inplace(double *x, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 1; i < w; i++)
		x[j*w+i] += x[j*w + i - 1];
	for (int i = 0; i < w; i++)
	for (int j = 1; j < h; j++)
		x[j*w+i] += x[(j-1)*w + i];
}

#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";
	int w, h, pd;
	double *x = iio_read_image_double_split(filename_in, &w, &h, &pd);
	for (int l = 0; l < pd; l++)
		compute_intimg_inplace(x + l*w*h, w, h);
	iio_write_image_double_split(filename_out, x, w, h, pd);
	return 0;
}
