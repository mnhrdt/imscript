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

#include "fail.c"
#include "xmalloc.c"


#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)



#ifndef USE_WISDOM
void evoke_wisdom(void) {}
void bequeath_wisdom(void) {}
#else//USE_WISDOM
#include "fftwisdom.c"
#endif//USE_WISDOM



// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2dfloat(fftwf_complex *fx, float *x, int w, int h)
{
	fftwf_complex *a = fftwf_malloc(w*h*sizeof*a);

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

//static void fft_2dfloatr2c(fftwf_complex *fx, float *x, int w, int h)
//{
//	//fftwf_complex *a = fftwf_malloc(w*h*sizeof*a);
//
//	//fprintf(stderr, "planning...\n");
//	evoke_wisdom();
//	fftwf_plan p = fftwf_plan_dft_r2c_2d(h, w, x, fx, FFTW_ESTIMATE);
//	bequeath_wisdom();
//	//fprintf(stderr, "...planned!\n");
//
//	//FORI(w*h) a[i] = x[i]; // complex assignment!
//	fftwf_execute(p);
//
//	fftwf_destroy_plan(p);
//	//fftwf_free(a);
//	fftwf_cleanup();
//}

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2dfloat(float *ifx,  fftwf_complex *fx, int w, int h)
{
	fftwf_complex *a = fftwf_malloc(w*h*sizeof*a);
	fftwf_complex *b = fftwf_malloc(w*h*sizeof*b);

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
		assert(cimagf(z) < 0.001);
	}
	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_free(b);
	fftwf_cleanup();
}


static void fft_direct(float *y, float *x, int w, int h, int pd)
{
	float *c = xmalloc(w*h*sizeof*c);
	fftwf_complex *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h)
			c[i] = x[i*pd + l];
		fft_2dfloat(gc, c, w, h);
		FORI(w*h) {
			y[2*(i*pd + l)+0] = crealf(gc[i]);
			y[2*(i*pd + l)+1] = cimagf(gc[i]);
		}
	}
	free(c);
	free(gc);
}

static void fft_inverse(float *y, float *x, int w, int h, int pd)
{
	int pdh = pd/2;
	assert(pd == 2*pdh);
	fftwf_complex *c = xmalloc(w*h*sizeof*c);
	float *gc = xmalloc(w*h*sizeof*gc);
	FORL(pdh) {
		FORI(w*h)
			c[i] = x[i*pd + 2*l] + I * x[i*pd+2*l+1];
		ifft_2dfloat(gc, c, w, h);
		FORI(w*h)
			y[i*pdh + l] = gc[i];
	}
	free(c);
	free(gc);
}

#ifndef OMIT_FFT_MAIN
int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s {1|-1} [in [out]]\n", *v);
		//                          0  1      2   3
		return EXIT_FAILURE;
	}
	int direction = atoi(v[1]);
	if (!direction) fail("bad direction");
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";

	int w, h, pd, pdout;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);

	if (direction < 0) {
		pdout = pd/2;
		if (2*pdout != pd) fail("can not ifft real data!");
	} else
		pdout = 2*pd;

	float *y = xmalloc(w*h*pdout*sizeof*y);

	if (direction < 0)
		fft_inverse(y, x, w, h, pd);
	else
		fft_direct(y, x, w, h, pd);

	iio_save_image_float_vec(out, y, w, h, pdout);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
#endif//OMIT_FFT_MAIN
