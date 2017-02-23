#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "iio.h"



#include "xmalloc.c"


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

	for (int i = 0; i < w*h; i++)
		a[i] = x[i]; // complex assignment!
	fftwf_execute(p);

	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_cleanup();
}

// if it finds any strange number, sets it to zero
static void normalize_float_array_inplace(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isnormal(x[i]))
			x[i] = 0;

}

static void dht(float *y, float *x, int w, int h)
{
	fftwf_complex *gc = xmalloc(w*h*sizeof*gc);
	fft_2dfloat(gc, x, w, h);
	for (int i = 0; i < w*h; i++)
		y[i] = (crealf(gc[i]) - cimagf(gc[i]))/sqrt(w*h);
	fftwf_free(gc);
}



#ifndef OMIT_DHT_MAIN
int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1      2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(in, &w, &h, &pd);
	normalize_float_array_inplace(x, w*h*pd);

	float *y = xmalloc(w*h*pd*sizeof*y);

	for (int i = 0; i < pd; i++)
		dht(y + i*w*h, x + i*w*h, w, h);

	iio_write_image_float_split(out, y, w, h, pd);
	free(x);
	free(y);
	return 0;
}
#endif//OMIT_DHT_MAIN
