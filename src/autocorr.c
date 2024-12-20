#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "iio.h"




#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

//#define π M_PI
//#define i∞ INFINITY

#include "fail.c"
#include "xmalloc.c"


#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)



#ifndef USE_WISDOM
static void evoke_wisdom(void) {}
static void bequeath_wisdom(void) {}
#else//USE_WISDOM
#include "fftwisdom.c"
#endif//USE_WISDOM



// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2dfloat(fftwf_complex *fx, float *x, int w, int h)
{
	fftwf_complex *a = fftwf_malloc(w*(long)h*sizeof*a);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_2d(h, w, a, fx,
						FFTW_FORWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = x[i]; // complex assignment!
	fftwf_execute(p);

	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_cleanup();
}


// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2dfloat(float *ifx,  fftwf_complex *fx, int w, int h)
{
	fftwf_complex *a = fftwf_malloc(w*(long)h*sizeof*a);
	fftwf_complex *b = fftwf_malloc(w*(long)h*sizeof*b);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_2d(h, w, a, b,
						FFTW_BACKWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = fx[i];
	fftwf_execute(p);
	float scale = 1.0/(w*h);
	FORI(w*h) {
		fftwf_complex z = b[i] * scale;
		ifx[i] = crealf(z);
		//assert(cimagf(z) < 0.001);
	}
	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_free(b);
	fftwf_cleanup();
}


// if it finds any strange number, sets it to zero
static void normalize_float_array_inplace(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isnormal(x[i]))
			x[i] = 0;

}



#define OMIT_PPSMOOTH_MAIN
#include "ppsmooth.c"


static int good_modulus(int n, int p)
{
	int r = n % p;
	r = r < 0 ? r + p : r;
	assert(r >= 0);
	assert(r < p);
	return r;
}


static void autocorr(
		float *y,  // output array to be filled-in
		float *xx, // input image
		int w,     // width
		int h,     // height
		float a,   // fractional autocorrelation coefficient
		float z,   // radius of central mask to zeroize
		int n,     // whether to normalize to zero mean
		int P      // whether to preprocess by periodic plus smooth
		)
{
	float *x = P ? xmalloc(w*h*sizeof*x) : xx;
	if (P) ppsmooth(x, xx, w, h);

	fftwf_complex *X = xmalloc(w*(long)h*sizeof*X);
	fft_2dfloat(X, x, w, h);
	if (a == 1)
		for (int i = 0; i < w*h; i++)
			X[i] = X[i] * conj(X[i]);
	else
		for (int i = 0; i < w*h; i++)
			X[i] = pow(cabs(X[i]), a);
	if (n) X[0] = 0;
	ifft_2dfloat(y, X, w, h);
	free(X);

	if (P) xfree(x);

	if (z > 0)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = good_modulus(i + w/2, w) - w/2;
		int jj = good_modulus(j + h/2, h) - h/2;
		if (hypot(ii, jj) < z)
			y[j*w + i] = 0;
	}
}

static
void autocorr_loc(float*y, float*x, int w,int h, float a, float z, int n, int L)
{
	// set background to 0
	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	// compute autocorr of all fitting square windows
	for (int oy = 0; oy+L < h; oy += L)
	for (int ox = 0; ox+L < w; ox += L)
	{
		float x_loc[L*L], y_loc[L*L], sy_loc[L*L];
		for (int j = 0; j < L; j++)
		for (int i = 0; i < L; i++)
			x_loc[j*L+i] = x[(j+oy)*w+i+ox];
		autocorr(y_loc, x_loc, L, L, a, z, n, 1);
		for (int j = 0; j < L; j++)
		for (int i = 0; i < L; i++)
		{
			int ii = (i + L/2) % L;
			int jj = (j + L/2) % L;
			sy_loc[j*L+i] = y_loc[jj*L+ii];
		}
		for (int j = 0; j < L; j++)
		for (int i = 0; i < L; i++)
			y[(j+oy)*w+i+ox] = sy_loc[j*L+i];
	}
}

static void autocorr_split(float *y,
		float *x, int w, int h, int pd, float a, float z, int n, int P)
{
	for (int i = 0; i < pd; i++)
		autocorr(y + w*h*i, x + w*h*i, w, h, a, z, n, P);
}

static void autocorr_loc_split(float *y,
		float *x, int w, int h, int pd, float a, float z, int n, int L)
{
	for (int i = 0; i < pd; i++)
		autocorr_loc(y + w*h*i, x + w*h*i, w, h, a, z, n, L);
}


typedef float (*extension_operator_float)(float*,int,int,int,int);

static float extend_float_image_periodic(float *x, int w, int h, int i, int j)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	return x[j*w+i];
}

static void fftshift(float *y, float *x, int w, int h)
{
	extension_operator_float p = extend_float_image_periodic;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[j*w+i] = p(x, w, h, i+w/2, j+h/2);

}

static void zeroize_inplace(float *x, int w, int h, float z)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int ii = good_modulus(i - w/2, w);
		int jj = good_modulus(j - h/2, h);
		if (hypot(ii, jj) < z)
			x[j*w + i] = 0;
	}
}

static void zeroize_inplace_split(float *x, int w, int h, int pd, float z)
{
	for (int i = 0; i < pd; i++)
		zeroize_inplace(x + w*h*i, w, h, z);
}

static void fftshift_split(float *y, float *x, int w, int h, int pd)
{
	for (int i = 0; i < pd; i++)
		fftshift(y + w*h*i, x + w*h*i, w, h);
}

#include "pickopt.c"
int main_autocorr(int c, char *v[])
{
	int   l = atoi(pick_option(&c, &v, "l", "0"));  // localization
	int   p = atoi(pick_option(&c, &v, "p", "0"));  // periodization
	int   s = atoi(pick_option(&c, &v, "s", "1"));  // shift 0 to center
	int   n = atoi(pick_option(&c, &v, "n", "1"));  // normalize average
	float a = atof(pick_option(&c, &v, "a", "1"));  // fractional alpha
	float R = atof(pick_option(&c, &v, "z", "0"));  // size of center mask
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = xmalloc(w*(long)h*pd*sizeof*y);
	normalize_float_array_inplace(x, w*h*pd);

	if (l)
		autocorr_loc_split(y, x, w, h, pd, a, R, n, l);
	else
		autocorr_split(y, x, w, h, pd, a, R, n, p);

	float *z = xmalloc(w*(long)h*pd*sizeof*z);
	if (s && !l)
		fftshift_split(z, y, w, h, pd);
	else
		z = y;


	iio_write_image_float_split(filename_out, z, w, h, pd);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_autocorr(c, v); }
#endif
