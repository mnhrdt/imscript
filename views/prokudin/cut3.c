// program to cut an image in three vertical parts

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

// x: input image
// w: width
// h: height
// r: output image to be filled-in with the upper part of image x
// g: output image to be filled-in with the middle part of image x
// b: output image to be filled-in with the bottom part of image x
void fill_three_parts(float *r, float *g, float *b, float *x, int w, int h)
{
	for (int j = 0; j < h/3; j++)
	for (int i = 0; i < w; i++)
	{
		r[j*w+i] = x[(j+0*h/3)*w+i];
		g[j*w+i] = x[(j+1*h/3)*w+i];
		b[j*w+i] = x[(j+2*h/3)*w+i];
	}
}

// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t%s in out1 out2 out3\n", *argv);
		//                          0 1  2    3    4
		return 1;
	}
	char *filename_in = argv[1];
	char *filename_out1 = argv[2];
	char *filename_out2 = argv[3];
	char *filename_out3 = argv[4];

	// read input image
	int w, h;
	float *in = iio_read_image_float(filename_in, &w, &h);

	// allocate space for the output images
	float *out1 = malloc(w*h*sizeof(float));
	float *out2 = malloc(w*h*sizeof(float));
	float *out3 = malloc(w*h*sizeof(float));

	// run the algorithm
	fill_three_parts(out1, out2, out3, in, w, h);

	// save the output images
	iio_save_image_float(filename_out1, out1, w, h/3);
	iio_save_image_float(filename_out2, out2, w, h/3);
	iio_save_image_float(filename_out3, out3, w, h/3);

	// cleanup and exit
	free(out1); free(out2); free(out3);
	free(in);
	return 0;
}
