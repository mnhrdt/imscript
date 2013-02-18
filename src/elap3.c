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



// evaluate the laplacian of image x at point i, j
inline static float laplacian(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = -4 * p(x, w, h, i  , j  )
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1);

	return r;
}


// returns the largest change performed all over the image
static float perform_one_iteration(float *x, int w, int h,
		int (*mask)[2], int nmask, float tstep)
{
	float maxupdate = 0;
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i;

		float new = x[idx] + tstep * laplacian(x, w, h, i, j);

		float update = fabs(x[idx] - new);
		if (update > maxupdate)
			maxupdate = update;

		x[idx] = new;
	}
	return maxupdate;
}


// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[2]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
		}
	assert(cx == nmask);

	*out_nmask = nmask;
	return mask;
}

// fill the holes of the image x using an harmonic function
void explicit_poisson_extension_with_initialization(
		float *y,        // output image
		float *x,        // input image (NAN values indicate holes)
		int w,           // image width
		int h,           // image height
		float timestep,  // time step for the numerical scheme
		int niter,       // number of iterations to run
		float *initialization
		)
{
	// build list of masked pixels
	int nmask, (*mask)[2] = build_mask(&nmask, x, w, h);

	// initialize the solution to the given data at the masked pixels
	for (int i = 0; i < w*h; i++)
		y[i] = isfinite(x[i]) ? x[i] : initialization[i];

	// do the requested iterations
	for (int i = 0; i < niter; i++)
	{
		float u = perform_one_iteration(y, w, h, mask, nmask, timestep);

		if (0 == i % 10)
			fprintf(stderr, "iter = %d, maxupdate = %g\n", i, u);
	}

	free(mask);
}


#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s TSTEP NITER data.png mask.png out.png ini.png\n", *argv);
		//0 1     2     3        4        5       6
		return 1;
	}
	float timestep = atof(argv[1]);
	int niter = atoi(argv[2]);
	char *filename_in = argv[3];
	char *filename_mask = argv[4];
	char *filename_out = argv[5];
	char *filename_ini = argv[6];

	int w[2], h[2];
	float *in = iio_read_image_float(filename_in, w, h);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	float *out = xmalloc(*w**h*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			in[i] = NAN;

	// If "filename_ini" is a number string, use this number as an
	// initialization constant.
	// If it is a real filename, use the corresponding image.
	float *init, constant = atof(filename_ini);
	if (0 != constant) {
		init = xmalloc(*w * *h * sizeof*init);
		for (int i = 0; i < *w**h; i++)
			init[i] = constant;
	} else {
		int ww, hh;
		init = iio_read_image_float(filename_ini, &ww, &hh);
		if (ww != *w || hh != *h)
			return fprintf(stderr, "initialization bad size");
	}

	explicit_poisson_extension_with_initialization(out,
					in, *w, *h, timestep, niter, init);

	iio_save_image_float(filename_out, out, *w, *h);

	return 0;
}
