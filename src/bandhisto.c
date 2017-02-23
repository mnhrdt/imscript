#include <stdio.h>
#include <math.h>

// input: an image of size wxh and b bands
// output: an image of size 256xb with the histogram of each band
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s img.tiff hist.tiff\n", *v);
		//                         0  1         2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	int nbins = 256;
	double H[pd*nbins];
	for (int i = 0; i < pd*nbins; i++)
		H[i] = 0;

	for (int l = 0; l < pd; l++)
	for (int i = 0; i < w*h; i++)
	{
		float v = x[w*h*l + i];
		if (!isfinite(v)) continue;
		int idx = floor(v);
		if (idx < 0) idx = 0;
		if (idx >= nbins) idx = nbins - 1;
		H[l*nbins+idx] += 1;
	}

	iio_write_image_double(filename_out, H, nbins, pd);

	return 0;
}
