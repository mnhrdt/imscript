// naive program to register two gray images
// method: find the translation that minimizes their L2 distance

#include <assert.h>
#include <math.h>
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
//
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
//
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


// zoom-out by a factor 2
//
// in: input image
// iw: input image width
// ih: input image height
// out: output image to be filled-in
// ow: output image width (supplied by the user)
// oh: output image height (supplied by the user)
//
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4];
		a[0] = getpixel_0(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel_0(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel_0(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel_0(in, iw, ih, 2*i+1, 2*j+1);
		out[ow*j + i] = (a[0] + a[1] + a[2] + a[3])/4;
	}
}


// evaluate the L2 distance between the image A and the image B translated by d
//
// A, B: input images
// w: width
// h: heigth
// d: translation vector
//
// return value: normalized L2 distance
//
float eval_displacement(float *A, float *B, int w, int h, int d[2])
{
	long double r = 0;
	int woff = w/16;
	int hoff = h/16;
	int npoints = (w-2*woff)*(h-2*hoff);
	for (int j = hoff; j < h-hoff; j++)
	for (int i = woff; i < w-woff; i++)
	{
		long double a = getpixel_0(A, w, h, i, j);
		long double b = getpixel_0(B, w, h, i-d[0], j-d[1]);
		r += (a - b) * (a - b)/(1.0*npoints);
	}
	r = sqrt(r);
	return r;
}


// recursive function to find the displacement vector that minimizes error
//
// A, B: input images
// w: width
// h: heigth
// scale: number of multi-scale recursions
// d: optimal displacement (output)
//
void find_displacement(int d[2], float *A, float *B, int w, int h, int scale)
{
	// find an initial rhough displacement d
	if (scale > 1) // call the function recursively
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *As = malloc(ws * hs * sizeof*As);
		float *Bs = malloc(ws * hs * sizeof*Bs);
		zoom_out_by_factor_two(As, ws, hs, A, w, h);
		zoom_out_by_factor_two(Bs, ws, hs, B, w, h);
		find_displacement(d, As, Bs, ws, hs, scale-1);
		free(As);
		free(Bs);
		d[0] *= 2;
		d[1] *= 2;
	} else {
		d[0] = 0;
		d[1] = 0;
	}

	// refine the rhough displacement by local optimization
	int bestn = -1, neig[9][2] = { {-1,-1},{-1,0},{-1,1},
		{0,-1},{0,0},{0,1}, {1,-1},{1,0},{1,1}, };
	float best = INFINITY;
	for (int n = 0; n < 9; n++)
	{
		int D[2] = {d[0] + neig[n][0], d[1] + neig[n][1]};
		float r = eval_displacement(A, B, w, h, D);
		if (r < best) {
			best = r;
			bestn = n;
		}
	}
	d[0] += neig[bestn][0];
	d[1] += neig[bestn][1];
	fprintf(stderr, "%dx%d: %d %d\n", w, h, d[0], d[1]);
}


// register two images
//
// w: width
// h: height
// left: left image
// right: right image
// out: right image after registration
//
void registration(float *out, float *left, float *right, int w, int h)
{
	int d[2];
	find_displacement(d, left, right, w, h, 10);
	apply_translation(out, d[0], d[1], right, w, h);
}


// main function
int main(int argc, char **argv)
{
	// process input arguments
	if (argc != 4) {
		fprintf(stderr, "usage:\n\t%s left right Tright\n", *argv);
		//                          0 1    2     3
		return 1;
	}
	char *filename_left = argv[1];
	char *filename_right = argv[2];
	char *filename_Tright = argv[3];

	// read input image
	int w, h;
	float *left = iio_read_image_float(filename_left, &w, &h);
	float *right = iio_read_image_float(filename_right, &w, &h);

	// allocate space for the output image
	float *out = malloc(w*h*sizeof(float));

	// run the algorithm
	registration(out, left, right, w, h);

	// save the output image
	iio_save_image_float(filename_Tright, out, w, h);

	// cleanup and exit
	free(left);
	free(right);
	free(out);
	return 0;
}
