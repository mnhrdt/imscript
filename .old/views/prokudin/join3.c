// merge three gray images into a RGB image
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

// w: width
// h: heigth
// r: input red channel
// g: input green channel
// b: input blue channel
// rgb: output color image, to be filled-in
void join_rgb(float *rgb, float *r, float *g, float *b, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		rgb[(j*w+i)*3+0] = r[j*w+i];
		rgb[(j*w+i)*3+1] = g[j*w+i];
		rgb[(j*w+i)*3+2] = b[j*w+i];
	}
}

// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t%s red green blue out\n", *argv);
		//                          0 1   2     3    4
		return 1;
	}
	char *filename_r = argv[1];
	char *filename_g = argv[2];
	char *filename_b = argv[3];
	char *filename_rgb = argv[4];

	// read input images
	int w[3], h[3];
	float *r = iio_read_image_float(filename_r, w+0, h+0);
	float *g = iio_read_image_float(filename_g, w+1, h+1);
	float *b = iio_read_image_float(filename_b, w+2, h+2);

	// check that the input images have the same size
	if (w[1]!=w[0] || w[2]!=w[0] || h[1]!=h[0] || h[2]!=h[0]) {
		fprintf(stderr, "input sizes mismatch(%d,%d)(%d,%d)(%d,%d)\n",
				w[0],h[0], w[1],h[1], w[2],h[2] );
		return 1;
	}

	// allocate space for the output image
	float *rgb = malloc(3*w[0]*h[0]*sizeof(float));

	// run the algorithm
	join_rgb(rgb, r, g, b, w[0], h[0]);

	// save the output image
	iio_save_image_float_vec(filename_rgb, rgb, w[0], h[0], 3);

	// cleanup and exit
	free(r); free(g); free(b);
	free(rgb);
	return 0;
}
