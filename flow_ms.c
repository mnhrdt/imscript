// tools for multi-scale optical flow {{{1
//
//
// There are two natural definitions for multi-scale optical flow:
//
// 1. Multi-scale by projection (can be used with any optical flow method)
// 2. Perturbative multi-scale (can be used with methods based on optimization)
//


// #includes {{{1
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "error.c"
#include "xmalloc.c"

// typedefs {{{1
typedef void (*generic_optical_flow)(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		int width, int height,
		void *data);

typedef void (*iterative_optical_flow)(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		float *in_u, float *in_v,
		int width, int height,
		void *data);

typedef float (*extension_operator_float)(
		float *image, int width, int height,
		int i, int j);

typedef float (*interpolation_operator_float)(
		float *image, int width, int height,
		float x, float y);

// utility functions {{{1
static float extend_float_image_constant(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j][i];
}

static float cell_interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float cell_interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	// return a;
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

static float cell_interpolate(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return cell_interpolate_nearest(a, b, c, d, x, y);
	//case 1: return marchi(a, b, c, d, x, y);
	case 2: return cell_interpolate_bilinear(a, b, c, d, x, y);
	default: return 0;
	}
	return -1;
}

static float interpolate_float_image_bilinearly(float *x, int w, int h,
		float i, float j)
{
	int ii = i;
	int jj = j;
	extension_operator_float p = extend_float_image_constant;
	float a = p(x, w, h, ii  , jj  );
	float b = p(x, w, h, ii  , jj+1);
	float c = p(x, w, h, ii+1, jj  );
	float d = p(x, w, h, ii+1, jj+1);
	return cell_interpolate(a, b, c, d, i-ii, j-jj, 2);
}

// image scaling {{{1
static void downscale_image(float *out, float *in,
		int outw, int outh, int inw, int inh,
		float scalestep)
{
	assert(scalestep == 2);
	fprintf(stderr, "downscale(%g): %dx%d => %dx%d\n", 
			scalestep, inw, inh, outw, outh);
	assert(2*outw >= inw);
	assert(2*outh >= inh);

	extension_operator_float p = extend_float_image_constant;

	for (int j = 0; j < outh; j++)
	for (int i = 0; i < outw; i++)
		out[outw*j+i] = 0.25 * ( p(in,inw,inh, 2*i, 2*j)
			       	+ p(in,inw,inh, 2*i+1, 2*j)
			       	+ p(in,inw,inh, 2*i, 2*j+1)
			       	+ p(in,inw,inh, 2*i+1, 2*j+1));
}

static void produce_downwards_pyramid(
		float **out_pyrx, int *out_pyrw, int *out_pyrh, int *out_pyrs,
		float *x, int w, int h,
		int nscales, float scalestep)
{
	assert(scalestep > 1);

	//// compute scale factors
	//float factor[nscales];
	//factor[0] = 1;
	//for (int i = 1; i < nscales; i++)
	//	factor[i] = factor[i-1] / scalestep;

	// compute pyramid sizes
	int pyrsize[nscales][3];
	pyrsize[0][0] = w;
	pyrsize[0][1] = h;
	for (int i = 1; i < nscales; i++) {
		pyrsize[i][0] = ceil(pyrsize[i-1][0] / scalestep);
		pyrsize[i][1] = ceil(pyrsize[i-1][1] / scalestep);
	}
	for (int i = 0; i < nscales; i++)
		pyrsize[i][2] = pyrsize[i][0] * pyrsize[i][1];

	// save requested outputt
	for (int i = 0; i < nscales; i++) {
		if (out_pyrw) out_pyrw[i] = pyrsize[i][0];
		if (out_pyrh) out_pyrh[i] = pyrsize[i][1];
		if (out_pyrs) out_pyrs[i] = pyrsize[i][2];
	}

	if (!out_pyrx) return;

	// alloc pyramid levels
	float *pyr[nscales];
	for (int i = 0; i < nscales; i++)
		pyr[i] = xmalloc(pyrsize[i][2] * sizeof(float));

	// fill initial level
	for (int i = 0; i < w*h; i++)
		pyr[0][i] = x[i];

	// propagate information to other levels
	if (out_pyrx)
	for (int i = 1; i < nscales; i++)
		downscale_image(pyr[i], pyr[i-1],
				pyrsize[i][0], pyrsize[i][1],
				pyrsize[i-1][0], pyrsize[i-1][1],
				scalestep);

	for (int i = 0; i < nscales; i++)
		if (out_pyrx) out_pyrx[i] = pyr[i];
}


static void upscale_flow(float *out, float *in,
		int outw, int outh, int inw, int inh,
		float scalestep)
{
	assert(scalestep == 2);
	assert(outw == 2*inw);
	assert(outh == 2*inh);

	float (*o)[outw] = (void*)out;
	float (*x)[inw] = (void*)in;

	// TODO: think about anchoring here, probably this is 0.5 pixels off!
	for (int j = 0; j < inh; j++)
	for (int i = 0; i < inw; i++)
	{
		o[2*j][2*i] = x[j][i];
		o[2*j+1][2*i] = (x[j][i] + x[j+1][i])/2;
		o[2*j][2*i+1] = (x[j][i] + x[j][i+1])/2;
		o[2*j+1][2*i+1] = (x[j][i]+x[j][i+1]+x[j+1][i]+x[j+1][i+1])/4;
	}

	for (int i = 0; i < outw * outh; i++)
		out[i] *= scalestep;
}

// multiscale flow {{{1
//void perturbative_multi_scale_optical_flow(float *u, float *v,
//		float *a, float *b, int w, int h,
//		iterative_optical_flow of, void *data,
//		int nscales, float scalestep)
//{

//
//	// compute scale factors
//	float factor[nscales];
//	factor[0] = 1;
//	for (int i = 1; i < nscales; i++)
//		factor[i] = factor[i-1] / scalestep;
//
//	// compute pyramid sizes
//	int pyrsize[nscales][2];
//	pyrsize[0][0] = w;
//	pyrsize[0][1] = h;
//	for (int i = 1; i < nscales; i++) {
//		pyrsize[i][0] = ceil(pyrsize[i-1][0] * factor[i]);
//		pyrsize[i][1] = ceil(pyrsize[i-1][1] * factor[i]);
//	}
//	for (int i = 0; i < nscales; i++)
//		pyrsize[i][2] = pyrsize[i][0] * pyrsize[i][1];
//
//	// build pyramid
//	float *pyr[nscales][4]; // (u, v, a, b)
//	for (int i = 0; i < nscales; i++)
//	for (int j = 0; j < 4; j++)
//		pyr[i][j] = xmalloc(pyrsize[i][2] * sizeof(float));
//	for (int i = 0; i < w*h; i++) {
//		pyr[0][2][i] = a[i];
//		pyr[0][3][i] = b[i];
//	}
//	for (int i = 1; i < nscales; i++) {
//		downscale_image(pyr[i][2], pyr[i-1][2],
//				pyrsize[i][0], pyrsize[i][1],
//				pyrsize[i-1][0], pyrsize[i-1][2], scalestep);
//		downscale_image(pyr[i][3], pyr[i-1][3],
//				pyrsize[i][0], pyrsize[i][1],
//				pyrsize[i-1][0], pyrsize[i-1][2], scalestep);
//	}
//
//
//	// run flow at top of pyramid, initialized by zero
//	float *tmpu = xmalloc(w * h * sizeof(float));
//	float *tmpv = xmalloc(w * h * sizeof(float));
//	int s = nscales-1;
//	for (int j = 0; j < pyrsize[s][2]; j++)
//		tmpu[j] = tmpv[j] = 0;
//	of(pyr[s][0], pyr[s][1], pyr[s][2], pyr[s][3], tmpu, tmpv, w, h, data);
//
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	s -= 1;
//	upscale_flow(pyr[s][2], pyr[s-1][2]);
//
//	// WORNG!!!! do start with the UPPER level, not the lower
//	// (XXX FIXME WRONG ERROR)
//
//	of(pyr[1][0], pyr[1][1], pyr[1][2], pyr[1][3], su, sv, w, h, data);
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	of(pyr[2][0], pyr[2][1], pyr[2][2], pyr[2][3], su, sv, w, h, data);
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	// ...
//}


//void generic_multi_scale_optical_flow(float *u, float *v,
//		float *a, float *b, int w, int h,
//		generic_optical_flow *of, void *data,
//		int nscales, float scalestep)
//{
//	assert(scalestep > 1);
//
//	// compute scale factors
//	float factor[nscales];
//	factor[0] = 1;
//	for (int i = 1; i < nscales; i++)
//		factor[i] = factor[i-1] / scalestep;
//
//	// compute pyramid sizes
//	int pyrsize[nscales][2];
//	pyrsize[0][0] = w;
//	pyrsize[0][1] = h;
//	for (int i = 1; i < nscales; i++) {
//		pyrsize[i][0] = ceil(pyrsize[i-1][0] * factor[i]);
//		pyrsize[i][1] = ceil(pyrsize[i-1][1] * factor[i]);
//	}
//	for (int i = 0; i < nscales; i++)
//		pyrsize[i][2] = pyrsize[i][0] * pyrsize[i][1];
//
//	// build pyramid
//	float *pyr[nscales][4]; // (u, v, a, b)
//	for (int i = 0; i < nscales; i++)
//	for (int j = 0; j < 4; j++)
//		pyr[i][j] = xmalloc(pyrsize[i][2] * sizeof(float));
//	for (int i = 0; i < w*h; i++) {
//		pyr[0][2][i] = a[i];
//		pyr[0][3][i] = b[i];
//	}
//	for (int i = 1; i < nscales; i++) {
//		downscale_image(pyr[i][2], pyr[i-1][2],
//				pyrsize[i][0], pyrsize[i][1],
//				pyrsize[i-1][0], pyrsize[i-1][2], scalestep);
//		downscale_image(pyr[i][3], pyr[i-1][3],
//				pyrsize[i][0], pyrsize[i][1],
//				pyrsize[i-1][0], pyrsize[i-1][2], scalestep);
//	}
//
//
//	// run flow at top of pyramid, initialized by zero
//	float *tmpu = xmalloc(w * h * sizeof(float));
//	float *tmpv = xmalloc(w * h * sizeof(float));
//	int s = nscales-1;
//	for (int j = 0; j < pyrsize[s][2]; j++)
//		tmpu[j] = tmpv[j] = 0;
//	of(pyr[s][0], pyr[s][1], pyr[s][2], pyr[s][3], tmpu, tmpv, w, h, data);
//
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	s -= 1;
//	upscale_flow(pyr[s][2], pyr[s-1][2]);
//
//	// WORNG!!!! do start with the UPPER level, not the lower
//	// (XXX FIXME WRONG ERROR)
//
//	of(pyr[1][0], pyr[1][1], pyr[1][2], pyr[1][3], su, sv, w, h, data);
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	of(pyr[2][0], pyr[2][1], pyr[2][2], pyr[2][3], su, sv, w, h, data);
//	// re-scale flow, and use it to initialize the next step of the pyramid
//	// ...
//}


#ifdef USE_MAIN

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s img step nscales outpat\n", *v);
		//                          0 1   2    3       4
		return EXIT_FAILURE;
	}
	char *filename_in = v[1];
	float scalestep = atof(v[2]);
	int nscales = atoi(v[3]);
	char *filepattern_out = v[4];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	float *pyrx[nscales];
	int pyrw[nscales], pyrh[nscales], pyrs[nscales];

	produce_downwards_pyramid(pyrx, pyrw, pyrh, pyrs,
			x, w, h, nscales, scalestep);

	for (int i = 0; i < nscales; i++) {
		char buf[0x100];
		snprintf(buf, 0x100, filepattern_out, i);
		printf("SCALE NUMBER %d: %dx%d (%d) => \"%s\"\n", i,
				pyrw[i], pyrh[i], pyrs[i], buf);
		iio_save_image_float(buf, pyrx[i], pyrw[i], pyrh[i]);
	}


	return EXIT_SUCCESS;
}

#endif//USE_MAIN

// vim:set foldmethod=marker:
