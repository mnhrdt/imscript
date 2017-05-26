// simpois: a simple Poisson solver

#include <math.h>
#include <stdlib.h>

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// evaluate the laplacian of image I at point i, j
static float laplacian(float *x, int w, int h, int i, int j)
{
	float r = -4 * getpixel_1(x, w, h, i  , j  )
		     + getpixel_1(x, w, h, i+1, j  )
		     + getpixel_1(x, w, h, i  , j+1)
		     + getpixel_1(x, w, h, i-1, j  )
		     + getpixel_1(x, w, h, i  , j-1);
	return r;
}

// perform one gauss-seidel iteration in-place on the image I (inside omega)
static void gauss_seidel_iteration(float *u, float *f, int w, int h,
		int (*omega)[2], int n_omega, float tstep)
{
	for (int p = 0; p < n_omega; p++)
	{
		int i = omega[p][0];
		int j = omega[p][1];
		int ij = j*w + i;

		float l = laplacian(u, w, h, i, j);
		float d = f ? f[ij] : 0;
		u[ij] = u[ij] + tstep * (l - d);
	}
}

// build a mask of the NAN positions on image "x" (region of interest)
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[2]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[2] = malloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
		}

	*out_nmask = nmask;
	return mask;
}

// iterative poisson solver with initialization
static void poisson_extension_with_init(
		float *u,        // output image
		float *f,        // interior Poisson data
		float *g,        // input image with boundary data (NAN = holes)
		int w,           // image width
		int h,           // image height
		float timestep,  // time step for the numerical scheme
		int niter,       // number of iterations to run
		float *initialization
		)
{
	// build the list of pixels inside the region of interest
	int n_omega, (*omega)[2] = build_mask(&n_omega, g, w, h);

	// initialize the solution to the given data outside the ROI
	for (int i = 0; i < w*h; i++)
		u[i] = isfinite(g[i]) ? g[i] : initialization[i];

	// perform the requested iterations
	for (int i = 0; i < niter; i++)
		gauss_seidel_iteration(u, f, w, h, omega, n_omega, timestep);

	// cleanup
	free(omega);
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_2(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	if (!out || !in) return;

	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = getpixel_1(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel_1(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel_1(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel_1(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_2(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = getpixel_1(in, iw, ih, i/2, j/2);
}

// recursive multi-scale poisson solver
void poisson_recursive(float *u, float *g, float *f, int w, int h,
		float tstep, int niter, int scale)
{
	float *init = malloc(w*h*sizeof*init);
	for (int i = 0 ; i < w*h; i++)
		init[i] = 0;

	if (scale > 1 && (w > 1 || h > 1)) {
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *gs = malloc(ws * hs * sizeof*gs);
		float *fs = f ? malloc(ws * hs * sizeof*fs) : f;
		float *us = malloc(ws * hs * sizeof*us);
		zoom_out_by_factor_2(gs, ws, hs, g, w, h);
		zoom_out_by_factor_2(fs, ws, hs, f, w, h);
		if (f) for (int i = 0; i < ws*hs; i++) fs[i] *= 3;
		poisson_recursive(us, gs, fs, ws, hs, tstep, niter, scale - 1);
		zoom_in_by_factor_2(init, w, h, us, ws, hs);
		free(gs); free(fs); free(us);
	}

	poisson_extension_with_init(u, f, g, w, h, tstep, niter, init);
	free(init);
}

// extension by Poisson equation of each channel of a color image
void poisson_solver_separable(float *out, float *in, float *dat,
		int w, int h, int pd, float tstep, int niter, int scale)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl  = in  + w*h*l;
		float *datl = dat ? dat + w*h*l : dat;
		poisson_recursive(outl, inl, datl, w, h, tstep, niter, scale);
	}
}


// comment the following line to disable the "main" function
//#define MAIN_IPOL_POISSON

#ifdef MAIN_IPOL_POISSON
#include <stdio.h>
#include "iio.h"      // library for image input/output
#include "pickopt.c"  // function for extracting named command line arguments

int main(int argc, char *argv[])
{
	// extract named arguments
	float tstep = atof(pick_option(&argc, &argv, "t", "0.45"));
	float niter = atof(pick_option(&argc, &argv, "n", "10"));
	float nscal = atof(pick_option(&argc, &argv, "s", "99"));
	char *filename_i = pick_option(&argc, &argv, "i", "-"); // stdin
	char *filename_o = pick_option(&argc, &argv, "o", "-"); // stdout
	char *filename_m = pick_option(&argc, &argv, "m", "");
	char *filename_f = pick_option(&argc, &argv, "f", "");

	// if any arguments are left, print a help message and quit
	if (argc > 1) {
		fprintf(stderr, "Usage:\n\t%s [options]\n", *argv);
		fprintf(stderr, "\nComputes approximate solutions of Poisson"
			" or Laplace equations for images\n");
		fprintf(stderr, "\nOptions with their default values:\n"
			"\t-i stdin   Input image with boundary data\n"
			"\t-f (zeros) Optional image with Poisson data term\n"
			"\t-m (zeros) Optional image with region of interest\n"
			"\t-o stdout  Output image\n"
			"\t-t 0.48    Time step for Gauss-Seidel iterations\n"
			"\t-n 10      Number of Gauss-Seidel iterations\n"
			"\t-s 99      Number of Multi-Scale octaves\n"
			"\nNote: NAN values in the input image"
			" are added to the region of interest\n");
		return 1;
	}

	// read input image ("boundary")
	int w, h, pd;
	float *img_i = iio_read_image_float_split(filename_i, &w, &h, &pd);

	// if requested, read data image
	float *img_f = NULL;
	if (*filename_f) {
		int w2, h2, pd2;
		img_f = iio_read_image_float_split(filename_f, &w2, &h2, &pd2);
		if (w != w2 || h != h2 || pd != pd2)
			return fprintf(stderr, "input sizes mismatch (-f)");
	}

	// if requested, read mask image
	float *img_m = NULL;
	if (*filename_m) {
		int w2, h2;
		img_m = iio_read_image_float(filename_m, &w2, &h2);
		if (w != w2 || h != h2)
			return fprintf(stderr, "input sizes mismatch (-m)");
	}

	// alloc space for output image
	float *out = malloc(w * h * pd * sizeof*out);

	// apply mask, if it exists
	if (img_m)
		for (int i = 0; i < w * h * pd; i++)
			if (img_m[i % (w*h)] > 0)
				img_i[i] = NAN;

	// run the algorithm
	poisson_solver_separable(out, img_i, img_f, w, h, pd,
			tstep, niter, nscal);

	// save the output image
	iio_save_image_float_split(filename_o, out, w, h, pd);

	// cleanup and exit
	free(out); free(img_i); free(img_m); free(img_f);
	return 0;
}
#endif//MAIN_IPOL_POISSON
