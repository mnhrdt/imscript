// program to quantize an image to 8-bits

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "iio.h"

// min: black point
// max: white point
// n: size of data
// x: input numbers
// y: output numbers, to be filled-in
void quantize(uint8_t *y, float *x, float min, float max, int n)
{
	for (int i = 0; i < n; i++)
	{
		float g = x[i];
		g = floor(255*(g - min)/(max - min));
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		y[i] = g;
	}
}

// main function
int main(int c, char *v[])
{
	// process input arguments
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s black white [in [out]]\n", *v);
		//                         0 1     2      3   4
		return 1;
	}
	float black = atof(v[1]);
	float white = atof(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	// read input image
	int w, h, pd;
	float *in = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	// allocate space for the output image
	uint8_t *out = malloc(w*h*pd);

	// run the algorithm
	quantize(out, in, black, white, w*h*pd);

	// save the output image
	iio_save_image_uint8_vec(filename_out, out, w, h, pd);

	// cleanup and exit
	free(out);
	free(in);
	return 0;
}
