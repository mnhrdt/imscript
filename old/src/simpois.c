#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


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

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel_1(float *I, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return I[i+j*w];
}

// evaluate the laplacian of image I at point i, j
static float laplacian(float *I, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = -4 * p(I, w, h, i  , j  )
		     + p(I, w, h, i+1, j  )
		     + p(I, w, h, i  , j+1)
		     + p(I, w, h, i-1, j  )
		     + p(I, w, h, i  , j-1);

	return r;
}

static float bilaplacian(float *x, int w, int h, int i, int j)
{
	float r = -4 * laplacian(x, w, h, i  , j  )
		     + laplacian(x, w, h, i+1, j  )
		     + laplacian(x, w, h, i  , j+1)
		     + laplacian(x, w, h, i-1, j  )
		     + laplacian(x, w, h, i  , j-1);

	return r;
}

static float fourlaplacian(float *x, int w, int h, int i, int j)
{
	float r = -4 * bilaplacian(x, w, h, i  , j  )
		     + bilaplacian(x, w, h, i+1, j  )
		     + bilaplacian(x, w, h, i  , j+1)
		     + bilaplacian(x, w, h, i-1, j  )
		     + bilaplacian(x, w, h, i  , j-1);

	return r;
}

// evaluate the laplacian of image I at point i, j
// (alternative function, compatible with neumann boundaries)
static float laplacian_neum(float *I, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float x = p(I, w, h, i  , j  );
	float a = p(I, w, h, i+1, j  );
	float b = p(I, w, h, i  , j+1);
	float c = p(I, w, h, i-1, j  );
	float d = p(I, w, h, i  , j-1);

	float r = 0;

	if (isfinite(a)) r += a -x;
	if (isfinite(b)) r += b -x;
	if (isfinite(c)) r += c -x;
	if (isfinite(d)) r += d -x;

	return r;
}

// perform one gauss-seidel iteration in-place on the data I
static void gauss_seidel_iteration(float *I, float *f, int w, int h,
		int (*omega)[2], int n_omega, float tstep)
{
	getpixel_operator op = laplacian_neum;
	if (tstep < 0)
	{
		op = bilaplacian;
		if (tstep < -1000)
		{
			op = fourlaplacian;
			tstep += 1000;
			tstep *= -1;
		}
	}

//#pragma omp parallel for
	for (int p = 0; p < n_omega; p++)
	{
		int i = omega[p][0];
		int j = omega[p][1];
		int ij = j*w + i;

		float l = op(I, w, h, i, j);
		float d = f ? f[ij] : 0;
		I[ij] = I[ij] + tstep * (l - d);
	}
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
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	if (!out || !in) return;
	getpixel_operator p = getpixel_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);

	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i, 2*j);
		a[1] = p(in, iw, ih, 2*i+1, 2*j);
		a[2] = p(in, iw, ih, 2*i, 2*j+1);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

// evaluate a bilinear cell at the given point
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = floor(p); // note: a cast to int fails when p<0
	int iq = floor(q);
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	if (!out || !in) return;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float x = (i - 0.5) / 2;
		float y = (j - 0.5) / 2;
		out[ow*j+i] = bilinear_interpolation(in, iw, ih, x, y);
		//out[ow*j+i] = getpixel_1(in, iw, ih, i/2, j/2);
	}
}

#include "cleant_cgpois.c"

#include "smapa.h"
SMART_PARAMETER(PMSFAC,3)
SMART_PARAMETER(PONLIT,0)

void poisson_rec(float *u, float *g, float *f, int w, int h,
		float tstep, int niter, int scale, int cgit)
{
	fprintf(stderr, "PREC %dx%d (niter,scale,cgit)=(%d %d %d)\n",
			w, h, niter, scale, cgit);
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1 && (w > 1 || h > 1))
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *gs = xmalloc(ws * hs * sizeof*gs);
		float *fs = f ? xmalloc(ws * hs * sizeof*fs) : f;
		float *us = xmalloc(ws * hs * sizeof*us);
		zoom_out_by_factor_two(gs, ws, hs, g, w, h);
		zoom_out_by_factor_two(fs, ws, hs, f, w, h);
		if (f) for (int i = 0; i < ws*hs; i++) fs[i] *= PMSFAC();
		poisson_rec(us, gs, fs, ws, hs, tstep, niter, scale-1, cgit);
		zoom_in_by_factor_two(init, w, h, us, ws, hs);
		free(gs);
		if(f) free(fs);
		free(us);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	if (PONLIT() && PONLIT() != w)
		niter = cgit = 0;

	poisson_extension_with_init(u, f, g, w, h, tstep, niter, init);
	free(init);

	if (cgit) { // if requested, refine by Conjugate Gradient
		float cg_eps = 1e-6;
		poisson_extension_by_cg(u, g, f, w, h, u, cgit, cg_eps);
	}
}

// extension by Poisson equation of each channel of a color image
void poisson_solver_separable(float *out, float *in, float *dat,
		int w, int h, int pd,
		float tstep, int niter, int scale, float cgrad)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl  = in  + w*h*l;
		float *datl = dat ? dat + w*h*l : dat;
		poisson_rec(outl, inl, datl, w, h, tstep, niter, scale, cgrad);
	}
}


#define MAIN_IPOL_POISSON
#ifdef MAIN_IPOL_POISSON
#include "iio.h"      // library for image input/output
#include "pickopt.c"  // function for extracting named command line arguments
int main_simpois(int argc, char *argv[])
{
	// extract named arguments
	float tstep = atof(pick_option(&argc, &argv, "t", "0.25"));
	float niter = atof(pick_option(&argc, &argv, "n", "10"));
	float nscal = atof(pick_option(&argc, &argv, "s", "99"));
	float cgrad = atof(pick_option(&argc, &argv, "c", "0"));
	char *filename_i = pick_option(&argc, &argv, "i", "-"); // stdin
	char *filename_o = pick_option(&argc, &argv, "o", "-"); // stdout
	char *filename_m = pick_option(&argc, &argv, "m", "");
	char *filename_f = pick_option(&argc, &argv, "f", "");
	if (!strcmp(filename_f, "-") && !strcmp(filename_i, "-"))
		filename_i = "NAN";

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
			"\t-t 0.25    Time step for Gauss-Seidel iterations\n"
			"\t-n 10      Number of Gauss-Seidel iterations\n"
			"\t-s 99      Number of Multi-Scale octaves\n"
			"\t-c 0       Number of Conjugate Gradient iterations\n"
		);
		fprintf(stderr, "\nNote: NAN values in the input image"
				" are added to the region of interest\n");
		return 1;
	}

	// read input image ("boundary")
	int w, h, pd;
	float *img_i = NULL;
	if (*filename_i && strcmp(filename_i, "NAN"))
		img_i = iio_read_image_float_split(filename_i, &w, &h, &pd);

	// if requested, read data image
	float *img_f = NULL;
	if (*filename_f) {
		int w2, h2, pd2;
		img_f = iio_read_image_float_split(filename_f, &w2, &h2, &pd2);
		if (!img_i) {
			w = w2; h = h2; pd = pd2;
			img_i = xmalloc(w * h * pd * sizeof*img_i);
			for (int i = 0; i < w*h*pd; i++) img_i[i] = NAN;
		}
		else if (w != w2 || h != h2 || pd != pd2)
			return fprintf(stderr, "input sizes mismatch (i,f)");
	}

	// if requested, read mask image
	float *img_m = NULL;
	if (*filename_m) {
		int w2, h2;
		img_m = iio_read_image_float(filename_m, &w2, &h2);
		if (w != w2 || h != h2)
			return fprintf(stderr, "input sizes mismatch (i,m)");
	}

	// alloc space for output image
	float *out = xmalloc(w * h * pd * sizeof*out);

	// apply mask, if it exists
	if (img_m)
		for (int i = 0; i < w * h * pd; i++)
			if (img_m[i % (w*h)] > 0)
				img_i[i] = NAN;

	// run the algorithm
	poisson_solver_separable(out, img_i, img_f, w, h, pd,
			tstep, niter, nscal, cgrad);

	// save the output image
	iio_write_image_float_split(filename_o, out, w, h, pd);

	// cleanup and exit
	free(out);
	free(img_i);
	if (img_m) free(img_m);
	if (img_f) free(img_f);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_simpois(c, v); }
#endif

#endif//MAIN_IPOL_POISSON
