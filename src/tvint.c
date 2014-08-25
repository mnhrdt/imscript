#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getsample_1(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l < 0) l = 0;
	if (l < pd) l = pd - 1;
	return x[(i+j*w)*pd+l];
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

// evaluate the divergence of a vector field F at point i, j
static float divergence(float *F, int w, int h, int i, int j)
{
	float x  = getsample_1(F, w, h, 2, i  , j  , 0);
	float y  = getsample_1(F, w, h, 2, i  , j  , 1);
	float xb = getsample_1(F, w, h, 2, i-1, j  , 0);
	float yb = getsample_1(F, w, h, 2, i  , j-1, 1);
	return x - xb + y - yb;
}

static void gradient(float *out_g, float *I, int w, int h, int i, int j)
{
	float I00   = getpixel_1(I, w, h, i  , j  );
	float I10 = getpixel_1(I, w, h, i+1, j  );
	float I01 = getpixel_1(I, w, h, i  , j+1);
	out_g[0] = I10 - I00;
	out_g[1] = I01 - I00;
}

//// perform one gauss-seidel iteration in-place on the data I
//static void gauss_seidel_iteration(float *I, int w, int h,
//		int (*omega)[2], int n_omega, float tstep)
//{
//	getpixel_operator op = laplacian_neum;
//	if (tstep < 0)
//		op = bilaplacian;
//
////#pragma omp parallel for
//	for (int p = 0; p < n_omega; p++)
//	{
//		int i = omega[p][0];
//		int j = omega[p][1];
//		int ij = j*w + i;
//
//		float l = op(I, w, h, i, j);
//		I[ij] = I[ij] + tstep * l;
//	}
//}

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

// perform one gauss-seidel iteration in-place on the data I
static void gauss_seidel_iteration(float *I, float *f, int w, int h,
		int (*omega)[2], int n_omega, float tstep)
{
	getpixel_operator op = laplacian_neum;
	if (tstep < 0)
		op = bilaplacian;

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
		if (isnan(g[i]))
			u[i] = initialization ? initialization[i] : 0;
		else
			u[i] = g[i];

	// perform the requested iterations
	for (int i = 0; i < niter; i++)
		gauss_seidel_iteration(u, f, w, h, omega, n_omega, timestep);

	// cleanup
	free(omega);
}


//// iterative minTV solver with initialization
//static void poisson_extension_with_init(
//		float *u,        // output image
//		float *g,        // input image with boundary data (NAN = holes)
//		int w,           // image width
//		int h,           // image height
//		float timestep,  // time step for the numerical scheme
//		int niter,       // number of iterations to run
//		float *initialization
//		)
//{
//	// build the list of pixels inside the region of interest
//	int n_omega, (*omega)[2] = build_mask(&n_omega, g, w, h);
//
//	// initialize the solution to the given data
//	for (int i = 0; i < w*h; i++)
//		if (isnan(g[i]))
//			u[i] = initialization[i]; // inside the ROI
//		else
//			u[i] = g[i];              // outside the ROI
//
//	// perform the requested iterations
//	for (int i = 0; i < niter; i++)
//		gauss_seidel_iteration(u, w, h, omega, n_omega, timestep);
//
//	// cleanup
//	free(omega);
//}

static void fill_p_from_initialization(float *p,
		float *g, float *init, int w, int h)
{
	float *v = xmalloc(w * h * sizeof*v);
	poisson_extension_with_init(v, init, g, w, h, 0.44, 100, NULL);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		gradient(p + 2 * (j*w+i), v, w, h, i, j);
	free(v);
}


// iterative minTV solver with initialization
static void tv_extension_with_init(
		float *u,        // output image
		float *g,        // input image with boundary data (NAN = holes)
		int w,           // image width
		int h,           // image height
		float tau,       // tau
		int niter,       // number of iterations to run
		float *initialization
		)
{
	// allocate space for the dual field and divergence field
	float *p = xmalloc(2 * w * h * sizeof*p);
	float *d = xmalloc(1 * w * h * sizeof*d);

	// initialize the dual field to zero everywhere (TODO: solve Poisson)
	fill_p_from_initialization(p, g, initialization, w, h);
	//for (int i = 0; i < 2*w*h; i++) p[i] = 0;

	// perform the requested iterations in the dual space
	for (int n = 0; n < niter; n++)
	{
		// fill the divergence field
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
			d[j*w+i] = divergence(p, w, h, i, j);

		// Chambolle projection
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			float G[2];
			gradient(G, d, w, h, i, j);
			float f = 1 + tau * hypot(g[0], g[1]);
			p[(j*w+i)*2 + 0] = ( p[(j*w+i)*2 + 0] + tau * G[0] ) / f;
			p[(j*w+i)*2 + 1] = ( p[(j*w+i)*2 + 1] + tau * G[1] ) / f;
		}
	}

	// recover the solution as the divergence of p
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		u[j*w+i] = divergence(p, w, h, i, j);

	// cleanup
	free(d);
	free(p);
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

#include "smapa.h"
SMART_PARAMETER(PONLIT,0)

void tvint_rec(float *u, float *g, int w, int h,
		float tstep, int niter, int scale)
{
	fprintf(stderr, "PREC %dx%d (niter,scale)=(%d %d)\n",
			w, h, niter, scale);
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1 && (w > 1 || h > 1))
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *gs = xmalloc(ws * hs * sizeof*gs);
		float *us = xmalloc(ws * hs * sizeof*us);
		zoom_out_by_factor_two(gs, ws, hs, g, w, h);
		tvint_rec(us, gs, ws, hs, tstep, niter, scale-1);
		zoom_in_by_factor_two(init, w, h, us, ws, hs);
		free(gs);
		free(us);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	if (PONLIT() && PONLIT() != w)
		niter = 0;

	tv_extension_with_init(u, g, w, h, tstep, niter, init);
	free(init);
}

// extension by TV minimization interpolation
void tv_interpolator_separable(float *out, float *in,
		int w, int h, int pd, int nscal, int niter, float tstep)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl  = in  + w*h*l;
		tvint_rec(outl, inl, w, h, tstep, niter, nscal);
	}
}


#define MAIN_IPOL_POISSON
#ifdef MAIN_IPOL_POISSON
#include "iio.h"      // library for image input/output
#include "pickopt.c"  // function for extracting named command line arguments
int main(int argc, char *argv[])
{
	// extract named arguments
	float tstep = atof(pick_option(&argc, &argv, "t", "0.25"));
	float niter = atof(pick_option(&argc, &argv, "n", "10"));
	float nscal = atof(pick_option(&argc, &argv, "s", "99"));
//	float cgrad = atof(pick_option(&argc, &argv, "c", "0"));
	char *filename_i = pick_option(&argc, &argv, "i", "-"); // stdin
	char *filename_o = pick_option(&argc, &argv, "o", "-"); // stdout
	char *filename_m = pick_option(&argc, &argv, "m", "");
//	char *filename_f = pick_option(&argc, &argv, "f", "");

	// if any arguments are left, print a help message and quit
	if (argc > 1) {
		fprintf(stderr, "Usage:\n\t%s [options]\n", *argv);
		fprintf(stderr, "\nComputes approximate solutions of Poisson"
			" or Laplace equations for images\n");
		fprintf(stderr, "\nOptions with their default values:\n"
			"\t-i stdin   Input image with boundary data\n"
//			"\t-f (zeros) Optional image with Poisson data term\n"
			"\t-m (zeros) Optional image with region of interest\n"
			"\t-o stdout  Output image\n"
//			"\t-t 0.25    Time step for Gauss-Seidel iterations\n"
//			"\t-n 10      Number of Gauss-Seidel iterations\n"
			"\t-s 99      Number of Multi-Scale octaves\n"
//			"\t-c 0       Number of Conjugate Gradient iterations\n"
		);
		fprintf(stderr, "\nNote: NAN values in the input image"
				" are added to the region of interest\n");
		return 1;
	}

	// read input image ("boundary")
	int w, h, pd;
	float *img_i = NULL;
	if (*filename_i)
		img_i = iio_read_image_float_split(filename_i, &w, &h, &pd);

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
	tv_interpolator_separable(out, img_i, w, h, pd, nscal, niter, tstep);

	// save the output image
	iio_save_image_float_split(filename_o, out, w, h, pd);

	// cleanup and exit
	free(out);
	free(img_i);
	free(img_m);
	return 0;
}
#endif//MAIN_IPOL_POISSON
