typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

// extrapolate by nearest value
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i + j*w];
}

void image_convolution_by_small_kernel(
		float *y,       // output image (to be filled-in)
		float *x,       // input image
		int w, int h,   // width and height of input and output images
		float *k,       // kernel
		int kw, int kh, // width and height of kernel
		int kp, int kq  // center coordinates of the kernel
		)
{
	getpixel_operator p = getpixel_0;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = 0;
		for (int jj = 0; jj < kh; jj++)
		for (int ii = 0; ii < kw; ii++)
		{
			int ci = i - kp + ii;
			int cj = j - kq + jj;
			a += p(k, kw, kh, ii, jj) * p(x, w, h, ci, cj);
		}
		y[j*w+i] = a;
	}
}


#ifdef CONVOLUTION_TEST_MAIN
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char **v)
{
	// data for available kernels
	float laplace[] = {3, 3, 1, 1, 0, 1, 0, 1, -4, 1, 0, 1, 0};
	float cx[] = {3, 1, 1, 0, -0.5, 0, 0.5};
	float cy[] = {1, 3, 0, 1, -0.5, 0, 0.5};
	float fx[] = {2, 1, 0, 0, -1, 1};
	float fy[] = {1, 2, 0, 0, -1, 1};
	float sobelx[] = {3, 3, 1, 1, -1, 0, 1, -2, 0, 2, -1, 0, 1};
	float sobely[] = {3, 3, 1, 1, -1, -2, -1, 0, 0, 0, 1, 2, 1};

	// process input arguments
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s filter [in [out]]\n", *v);
		//                          0 1       2   3
		return 1;
	}
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";
	float *filter = NULL;
	if (0 == strcmp(v[1], "laplace")) filter = laplace;
	if (0 == strcmp(v[1], "cx"     )) filter = cx;
	if (0 == strcmp(v[1], "cy"     )) filter = cy;
	if (0 == strcmp(v[1], "fx"     )) filter = fx;
	if (0 == strcmp(v[1], "fy"     )) filter = fy;
	if (0 == strcmp(v[1], "sobelx" )) filter = sobelx;
	if (0 == strcmp(v[1], "sobely" )) filter = sobely;
	if (!filter) {
		fprintf(stderr, "filters = cx, fx, sobelx, laplace, cy, ...\n");
		return 1;
	}

	// fill kernel parameters
	int k_width   = filter[0];
	int k_height  = filter[1];
	int k_centerx = filter[2];
	int k_centery = filter[3];
	float *k_data = filter + 4;

	// prepare input and output images
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *y = malloc(w*h*sizeof*y);

	// compute
	image_convolution_by_small_kernel(y, x, w, h, k_data,
			k_width, k_height, k_centerx, k_centery);

	// save result
	iio_save_image_float(filename_out, y, w, h);

	// cleanup
	free(x);
	free(y);
	return 0;
}
#endif//CONVOLUTION_TEST_MAIN
