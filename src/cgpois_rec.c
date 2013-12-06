#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


#include "conjugate_gradient.c"

// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[3]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[3] = xmalloc(w*h*3*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			mask[cx][2] = j*w+i;
			cx += 1;
		}
	assert(cx == nmask);

	*out_nmask = nmask;
	return mask;
}

static void invert_mask(int *im, int (*mask)[3], int nmask, int w, int h)
{
	for (int i = 0; i < w*h; i++)
		im[i] = -1;
	for (int i = 0; i < nmask; i++)
	{
		int idx = mask[i][2];// + w*mask[i][1];
		assert(idx >= 0);
		assert(idx < w*h);
		im[idx] = i;
	}
}

struct cgpois_state {
	int w, h, (*mask)[3], nmask, *invmask;
	float *boundary_data;
	float *interior_data;
};

typedef float (*fancy_getpixel_operator)(double*x,void*,int,int);

static float masked_getpixel_1(double *x, void *ee, int i, int j)
{
	struct cgpois_state *e = ee;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= e->w) i = e->w-1;
	if (j >= e->h) j = e->h-1;
	int idx = j*e->w + i;
	assert(idx >= 0);
	assert(idx < e->w * e->h);
	//float *y = e->invmask[idx] < 0 ? e->original_image : x;
	//return y[idx];
	int pidx = e->invmask[idx];
	if (pidx < 0) {
		//assert(isfinite(e->original_image[idx]));
		return e->boundary_data[idx];
	} else {
		assert(pidx < e->nmask);
		//assert(isnan(e->original_image[idx]));
		return x[pidx];
	}
}

static double evaluate_laplacian_at(double *x, int p, void *ee)
{
	struct cgpois_state *e = ee;
	assert(p >= 0);
	assert(p < e->nmask);
	int i = e->mask[p][0];
	int j = e->mask[p][1];

	fancy_getpixel_operator g = masked_getpixel_1;

	double r = -4 * g(x, ee, i  , j  )
		     + g(x, ee, i+1, j  )
		     + g(x, ee, i  , j+1)
		     + g(x, ee, i-1, j  )
		     + g(x, ee, i  , j-1);
	return r;
}

static void minus_operator(double *y, double *x, int n, void *ee)
{
	for (int p = 0; p < n; p++)
		y[p] = -evaluate_laplacian_at(x, p, ee);
}

#include "smapa.h"
SMART_PARAMETER(CG_MAXIT,-1)
SMART_PARAMETER(CG_EPS,-1)


void poisson_solver_with_init(float *out, float *in, float *dat, int w, int h,
		float *init)
{
	// build list of masked pixels
	int nmask, (*mask)[3] = build_mask(&nmask, in, w, h);
	if (!nmask) {
		for (int i = 0; i < w*h; i++)
			out[i] = isfinite(in[i]) ? in[i] : init[i];
		return;
	}
	int *invmask = xmalloc(w*h*sizeof(int));
	invert_mask(invmask, mask, nmask, w, h);

	// define the linear map A=laplacian_operator
	struct cgpois_state e[1];
	e->w = w;
	e->h = h;
	e->mask = mask;
	e->invmask = invmask;
	e->nmask = nmask;
	e->boundary_data = xmalloc(w*h*sizeof(double));
	e->interior_data = xmalloc(w*h*sizeof(double));
	for (int i = 0; i < w*h; i++)
		e->boundary_data[i] = in[i];
	for (int i = 0; i < w*h; i++)
		e->interior_data[i] = 0;//dat[i];

	// fill-in the independent term b
	double *b = xmalloc(nmask * sizeof(double));
	double *tmp = xmalloc(nmask * sizeof(double));
	for (int p = 0; p < nmask; p++)
		tmp[p] = 0;
	minus_operator(b, tmp, nmask, e);
	free(tmp);
	for (int p = 0; p < nmask; p++)
		b[p] = -dat[mask[p][0]+w*mask[p][1]] - b[p];
	for (int i = 0; i < w*h; i++)
		if (e->invmask[i] < 0)
			e->boundary_data[i] = 0;

	// compute the solution
	double *solution = xmalloc(nmask * sizeof(double));
	double *initialization = xmalloc(nmask * sizeof(double));
	for (int i = 0; i < nmask; i++)
		initialization[i] = init[mask[i][0]+w*mask[i][1]];
	//conjugate_gradient(solution, minus_laplacian_operator, b, nmask, e);
	int cg_maxit = CG_MAXIT() >= 0 ? CG_MAXIT() : nmask;
	float cg_eps = CG_EPS() >= 0 ? CG_EPS() : 1e-6;
	fancy_conjugate_gradient(solution, minus_operator, b, nmask,
					e, initialization, cg_maxit, cg_eps);

	// copy the solution to its place
	for (int i = 0; i < w*h; i++)
		out[i] = in[i];
	for (int p = 0; p < nmask; p++) {
		if (nmask < 33)
			fprintf(stderr, "sol[%d] = %g\n", p, solution[p]);
		out[mask[p][2]] = solution[p];
	}

	free(e->interior_data);
	free(e->boundary_data);
	free(solution);
	free(mask);
	free(invmask);
	free(b);
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


// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
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

static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
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
	getpixel_operator p = getpixel_1;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float x = (i - 0.5)/2;
		float y = (j - 0.5)/2;
		out[ow*j+i] = bilinear_interpolation(in, iw, ih, x, y);
		//out[ow*j+i] = p(in, iw, ih, round(x), round(y));
	}
}

void poisson_recursive(float *out, float *in, float *dat, int w, int h,
		int scale)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins  = xmalloc(ws * hs * sizeof*ins);
		float *dats = xmalloc(ws * hs * sizeof*dats);
		float *outs = xmalloc(ws * hs * sizeof*outs);
		zoom_out_by_factor_two(ins, ws, hs, in, w, h);
		zoom_out_by_factor_two(dats, ws, hs, dat, w, h);
		for (int i = 0; i < ws*hs; i++)
			dats[i] *= 4;
		poisson_recursive(outs, ins, dats, ws, hs, scale-1);
		zoom_in_by_factor_two(init, w, h, outs, ws, hs);

		free(ins);
		free(dats);
		free(outs);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	poisson_solver_with_init(out, in, dat, w, h, init);
	free(init);
}


#include "iio.h"

SMART_PARAMETER(NSCALES,1)

int main(int argc, char *argv[])
{
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s boundary.png data.png mask.png out.png\n", *argv);
		//        0 1            2        3        4
		return 1;
	}
	char *filename_inpu = argv[1];
	char *filename_data = argv[2];
	char *filename_mask = argv[3];
	char *filename_out = argv[4];

	int w[3], h[3];
	float *inpu = iio_read_image_float(filename_inpu, w, h);
	float *data = iio_read_image_float(filename_data, w+1, h+1);
	float *mask = iio_read_image_float(filename_mask, w+2, h+2);
	if (w[0] != w[1] || h[0] != h[1] || w[0] != w[2] || h[0] != h[2])
		return fprintf(stderr, "input image files sizes mismatch");
	float *out = xmalloc(*w**h*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			inpu[i] = NAN;

	int nscales = NSCALES();
	poisson_recursive(out, inpu, data, *w, *h, nscales);

	iio_save_image_float(filename_out, out, *w, *h);

	return 0;
}
