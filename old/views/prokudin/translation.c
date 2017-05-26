// translate an image by an integral number of pixels

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

// auxiliary function to get the value of an image at any point (i,j)
// (points outisde the original domain get the value 0)
//
// x: image data
// w: width
// h: height
// i: horizontal position
// j: vertical position
//
// return value: color of the requested pixel
float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	else
		return x[j*w+i];
}


// apply a translation to the given image
//
// in: input image
// w: width
// h: height
// dx: horizontal displacement
// dy: vertical displacement
// out: output image, to be filled-in
void apply_translation(float *out, int dx, int dy, float *in, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = i - dx;
		int jj = j - dy;
		out[j*w+i] = getpixel_0(in, w, h, ii, jj);
	}
}


// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t%s dx dy in out\n", *argv);
		//                          0 1   2 3  4
		return 1;
	}
	int dx = atoi(argv[1]);
	int dy = atoi(argv[2]);
	char *filename_in = argv[3];
	char *filename_out = argv[4];

	// read input image
	int w, h;
	float *in = iio_read_image_float(filename_in, &w, &h);

	// allocate space for output image
	float *out = malloc(w*h*sizeof(float));

	// run the algorithm
	apply_translation(out, dx, dy, in, w, h);

	// save output image
	iio_save_image_float(filename_out, out, w, h);

	// cleanup and exit
	free(in);
	free(out);
	return 0;
}
