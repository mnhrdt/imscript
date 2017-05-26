#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// utility function that always returns a valid pointer to memory
static void *xmalloc(size_t n)
{
	void *new = malloc(n);
	if (!new)
	{
		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
		exit(1);
	}
	return new;
}


// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
inline static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w];
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
inline static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

inline static float ricciop(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = -4 * log(p(x, w, h, i  , j  ))
		     + log(p(x, w, h, i+1, j  ))
		     + log(p(x, w, h, i  , j+1))
		     + log(p(x, w, h, i-1, j  ))
		     + log(p(x, w, h, i  , j-1));

	return r/2;
}

inline static float ricciop2(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float gx = p(x, w, h, i+1, j) - p(x, w, h, i-1, j);
	float gy = p(x, w, h, i, j+1) - p(x, w, h, i, j-1);
	float n = hypot(gx, gy);
	n = hypot(0.01, n);

	float l = -4 * p(x, w, h, i  , j  )
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1);

	return l*n/255;
}


// returns the largest change performed all over the image
static float perform_one_iteration(float *y, float *x, int w, int h,
		float tstep)
{
	float maxupdate = 0;

	for (int j = 0; j < w; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;

		float new = x[idx] + tstep * ricciop2(x, w, h, i, j);

		float update = fabs(x[idx] - new);
		if (update > maxupdate)
			maxupdate = update;

		y[idx] = new;
	}
	return maxupdate;
}



void isoricci(
		float *y,        // output image
		float *x,        // input image
		int w,           // image width
		int h,           // image height
		float timestep,  // time step for the numerical scheme
		int niter        // number of iterations to run
		)
{
	// do the requested iterations
	for (int i = 0; i < niter; i++)
	{
		float u = perform_one_iteration(y, x, w, h, timestep);
		for (int j = 0; j < w*h; j++)
			x[j] = y[j];

		if (0 == i % 10)
			fprintf(stderr, "iter = %d, maxupdate = %g\n", i, u);
	}
}


#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s TSTEP NITER in.png out.png\n", *argv);
		//        0 1     2     3      4
		return 1;
	}
	float timestep = atof(argv[1]);
	int niter = atoi(argv[2]);
	char *filename_in = argv[3];
	char *filename_out = argv[4];

	int w[2], h[2];
	float *in = iio_read_image_float(filename_in, w, h);
	float *out = xmalloc(*w**h*sizeof*out);

	isoricci(out, in, *w, *h, timestep, niter);

	iio_write_image_float(filename_out, out, *w, *h);

	return 0;
}
