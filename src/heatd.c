float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

// anisotropic laplacian
float laplacian(
		float *x,        // image
		float *A,        // diffusion matrix field (a,b,c)
		int w,           // image width
		int h,           // image height
		int i,           // i-position of the point
		int j            // j-position of the pointj
		)
{
	// default diffusion matrix
	float a = 1;
	float b = 0;
	float c = 1;

	// if provided, use another matrix
	if (A)
	{
		a = A[3*(j*w+i)+0];
		b = A[3*(j*w+i)+1];
		c = A[3*(j*w+i)+2];
	}

	// build scheme
	float s[9] = {
		b/2  , c        , -b/2 ,
		a  , -2*a-2*c , a  ,
		-b/2 , c        , b/2
	};

	// extract nine values
	float v[9] = {
		getpixel(x, w, h, i-1, j-1),
		getpixel(x, w, h, i  , j-1),
		getpixel(x, w, h, i+1, j-1),
		getpixel(x, w, h, i-1, j  ),
		getpixel(x, w, h, i  , j  ),
		getpixel(x, w, h, i+1, j  ),
		getpixel(x, w, h, i-1, j+1),
		getpixel(x, w, h, i  , j+1),
		getpixel(x, w, h, i+1, j+1),
	};

	// apply scheme
	float r = 0;
	for (int k = 0; k < 9; k++)
		r += s[k] * v[k];
	return r/4;
}

void linear_diffusion_step(
		float *y,        // output image
		float *x,        // input image
		int w,           // image width
		int h,           // image height
		float *A,        // diffusion matrix field (a,b,c)
		float t          // timestep
		)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[j*w+i] = x[j*w+i] + t * laplacian(x, A, w, h, i, j);
}

void linear_diffusion(
		float *y,        // output image
		float *x,        // input image
		int w,           // image width
		int h,           // image height
		float *A,        // diffusion matrix field (a,b,c)
		int n,           // number of diffusion steps
		float t          // timestep
		)
{
	for (int i = 0; i < n; i++)
	{
		linear_diffusion_step(y, x, w, h, A, t);
		for (int j = 0; j < w*h; j++)
			x[j] = y[j];
	}
}

// note (fills "x" with the n-1th iteration)
void linear_diffusion_separable(
		float *y,        // output image
		float *x,        // input image
		int w,           // image width
		int h,           // image height
		int pd,          // number of channels (pixel dimension)
		float *A,        // diffusion matrix field (a,b,c)
		int n,           // number of diffusion steps
		float t          // timestep
		)
{
	for (int i = 0; i < pd; i++)
		linear_diffusion(y + i*w*h, x + i*w*h, w, h, A, n, t);
}


#include <stdio.h> // fprintf
#include <stdlib.h> // atoi, atof
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	// extract named options
	int  nsteps = atoi(pick_option(&c, &v, "n", "1"));
	float tstep = atof(pick_option(&c, &v, "t", "0.25"));
	char *filename_g = pick_option(&c, &v, "g", "");
	char *filename_A = pick_option(&c, &v, "A", "");
	char *filename_m = pick_option(&c, &v, "m", "");
	char *matrix_D   = pick_option(&c, &v, "D", "");

	// extract positional arguments
	if (c != 1 && c != 2 && c != 3)
		return fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
	//                                         0  1   2
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";


	// read input images
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *A = NULL;
	if (*filename_A) {
		int ww, hh, ppdd;
		A = iio_read_image_float_vec(filename_A, &ww, &hh, &ppdd);
		if (ww != w || hh != h || ppdd != 3)
			return fprintf(stderr, "ERROR: bad A dimensions"
				"(%d %d %d) ((%d %d))\n", ww, hh, ppdd, w, h);
	}
	float *y = malloc(w * h * pd * sizeof*y);

	// call algorithm
	linear_diffusion_separable(y, x, w, h, pd, A, nsteps, tstep);

	// save result and quit
	iio_write_image_float_split(filename_out, y, w, h, pd);
	return 0;
}
