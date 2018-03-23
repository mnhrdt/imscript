#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "smapa.h"
SMART_PARAMETER(NSCALES,20)
SMART_PARAMETER(TSTEP,0.25)
SMART_PARAMETER(NITER,10)

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = getpixel(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				//m += a[k];
				//cx += 1;
				m = a[k];
				cx = 1;
				break;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two_neum(float *out, int ow, int oh,
		float *in, int iw, int ih, float *neumask)
	// neum has size iwxih
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = getpixel(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				//m += a[k];
				//cx += 1;
				m = a[k];
				cx = 1;
				break;
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
	int ip = p;
	int iq = q;
	float a = getpixel(x, w, h, ip  , iq  );
	float b = getpixel(x, w, h, ip+1, iq  );
	float c = getpixel(x, w, h, ip  , iq+1);
	float d = getpixel(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = bilinear_interpolation(in, iw, ih,
						(i-0.5)/2, (j-0.5)/2);
}

static void zoom_in_by_factor_two_neum(float *out, int ow, int oh,
		float *in, int iw, int ih, float *neum)
	// neum has size iwxih
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = bilinear_interpolation(in, iw, ih,
						(i-0.5)/2, (j-0.5)/2);
}

static float laplacian_neum(float *I, float *n, int w, int h, int i, int j)
{
	// image values
	float x = getpixel(I, w, h, i  , j  );
	float a = getpixel(I, w, h, i+1, j  );
	float b = getpixel(I, w, h, i  , j+1);
	float c = getpixel(I, w, h, i-1, j  );
	float d = getpixel(I, w, h, i  , j-1);

	// neumask values
	int nx = 0 < getpixel(n, w, h, i  , j  );
	int na = 0 < getpixel(n, w, h, i+1, j  );
	int nb = 0 < getpixel(n, w, h, i  , j+1);
	int nc = 0 < getpixel(n, w, h, i-1, j  );
	int nd = 0 < getpixel(n, w, h, i  , j-1);

	// computation
	float r = 0;
	if (nx)
	{
		if (!na && isfinite(a)) r += a-x;
		if (!nb && isfinite(b)) r += b-x;
		if (!nc && isfinite(c)) r += c-x;
		if (!nd && isfinite(d)) r += d-x;
	}
	return r;
}

// perform one gauss-seidel iteration in-place on the data I
static void gauss_seidel_neumann(
		float *I,
		float *n,
		int w,
		int h,
		int (*omega)[2],
		int n_omega,
		float tstep
		)
{
//#pragma omp parallel for
	for (int p = 0; p < n_omega; p++)
	{
		int i = omega[p][0];
		int j = omega[p][1];
		int ij = j*w + i;

		float l = laplacian_neum(I, n, w, h, i, j);
		if (fabs(l) > 0.2) fprintf(stderr, "\tij=(%d %d){%g} l=%g\n",i,j,I[ij],l);
		float d = 0;//f ? f[ij] : 0;
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
	int (*mask)[2] = malloc(w*h*2*sizeof(int)), cx = 0;
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

static void run_a_few_poisson_iterations(
		float *u,        // output image
		float *g,        // input image with boundary data (NAN = holes)
		int w,           // image width
		int h,           // image height
		float *n,        // neumann boundary mask
		float timestep,  // time step for numerical scheme
		int niter,       // number of iterations to run
		float *initialization
		)
{
	fprintf(stderr, "r_p_i %d %d {%d}\n", w, h, niter);
	// build the list of pixels inside the region of interest
	int n_omega, (*omega)[2] = build_mask(&n_omega, g, w, h);

	// initialize the solution to the given data outside the ROI
	for (int i = 0; i < w*h; i++)
		u[i] = isfinite(g[i]) ? g[i] : initialization[i];

	// perform the requested iterations
	for (int i = 0; i < niter; i++)
	{
		fprintf(stderr, "\titer %d {nomega=%d}\n", i, n_omega);
		gauss_seidel_neumann(u, n, w, h, omega, n_omega, timestep);
	}

	// cleanup
	free(omega);
}

void simplest_inpainting_neum(float *out, float *in, float *neum,
		int w, int h, int scale)
{
	fprintf(stderr, "s_i_n %d %d\n", w, h);
	float *init = malloc(w * h * sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins  = malloc(ws * hs * sizeof*ins);
		float *outs = malloc(ws * hs * sizeof*outs);
		float *neus  = malloc(ws * hs * sizeof*neus);
		zoom_out_by_factor_two(neus, ws, hs, neum, w, h);
		zoom_out_by_factor_two_neum(ins, ws, hs, in, w, h, neum);
		simplest_inpainting_neum(outs, ins, neus, ws, hs, scale - 1);
		zoom_in_by_factor_two_neum(init, w, h, outs, ws, hs, neum);
		free(ins);
		free(outs);
		free(neus);
	}
	else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 33;
	}
#if 0
	for (int i = 0 ; i < w*h; i++)
		out[i] = isfinite(in[i]) ? in[i] : init[i];
#else
	float tstep = TSTEP();
	int niter = NITER();
	run_a_few_poisson_iterations(out, in, w, h, neum, tstep, niter, init);
#endif
	free(init);
}

void simplest_inpainting_separable_neum(float *out, float *in, float *neum,
		int w, int h, int pd, int nscales)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		simplest_inpainting_neum(outl, inl, neum, w, h, nscales);
	}
}


#include "iio.h"
int main(int argc, char *argv[])
{
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t"
		"%s data.png mask.png neumann.png out.png\n", *argv);
		//0 1        2        3           4
		return 1;
	}
	char *filename_in = argv[1];
	char *filename_mask = argv[2];
	char *filename_neum = argv[3];
	char *filename_out = argv[4];

	int nscales = NSCALES();

	int w[3], h[3], pd;
	float *in = iio_read_image_float_split(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	float *neum = iio_read_image_float(filename_mask, w+2, h+2);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	if (w[0] != w[2] || h[0] != h[2])
		return fprintf(stderr, "image and neum file size mismatch");
	float *out = malloc(*w * *h *pd * sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			for (int l = 0; l < pd; l++)
				in[*w**h*l+i] = NAN;

	simplest_inpainting_separable_neum(out, in, neum, *w, *h, pd, nscales);

	iio_write_image_float_split(filename_out, out, *w, *h, pd);

	return 0;
}
