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

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)


static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "out of memory\n"));
	return new;
}


#ifndef USE_WISDOM
static void evoke_wisdom(void) {}
static void bequeath_wisdom(void) {}
#else//USE_WISDOM
#include "fftwisdom.c"
#endif//USE_WISDOM


static void dct_2dfloat(float *fx, float *x, int w, int h)
{
	float normalization_factor = sqrt(4*(w-1)*(h-1));
	float *a = fftwf_malloc(w*h*sizeof*a);
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_r2r_2d(h, w, a, fx,
			FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
	bequeath_wisdom();
	FORI(w*h) a[i] = x[i] / normalization_factor;
	fftwf_execute(p);
	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_cleanup();
}

static void dct(float *y, float *x, int w, int h, int pd)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h)
			c[i] = x[i*pd + l];
		dct_2dfloat(gc, c, w, h);
		FORI(w*h)
			y[i*pd + l] = gc[i];
	}
	free(c);
	free(gc);
}

// if it finds any strange number, sets it to zero
static void normalize_float_array_inplace(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isnormal(x[i]))
			x[i] = 0;

}

int main_dct(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	normalize_float_array_inplace(x, w*h*pd);

	float *y = xmalloc(w*h*pd*sizeof*y);

	dct(y, x, w, h, pd);

	iio_write_image_float_vec(out, y, w, h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_dct(c, v); }
#endif
