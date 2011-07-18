// gaussian blur of a 2D image

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

// this function prints an error message and aborts the program
static void error(const char *fmt, ...)
{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#ifdef NDEBUG
	exit(-1);
#else
	exit(*(int *)0x43);
#endif
}


// this function always returns a valid pointer
static void *xmalloc(size_t size)
{
	if (size == 0)
		error("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
	{
		double sm = size / (0x100000 * 1.0);
		error("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return new;
}


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

//static void pointwise_complex_rmultiplication(fftwf_complex *w,
//		fftwf_complex *z, float *x, int n)
//{
//	FORI(n)
//		w[i] = z[i] * x[i];
//}

static void pointwise_complex_multiplication(fftwf_complex *w,
		fftwf_complex *z, fftwf_complex *x, int n)
{
	FORI(n)
		w[i] = z[i] * x[i];
}

static void fill_2d_gaussian_image(float *gg, int w, int h, float inv_s)
{
	float (*g)[w] = (void *)gg;
	float alpha = inv_s*inv_s/(M_PI);

	FORJ(h) FORI(w) {
		float x = i < w/2 ? i : i - w;
		float y = j < h/2 ? j : j - h;
		float r = hypot(x, y);
		g[j][i] = alpha * exp(-r*r*inv_s*inv_s);
	}

	// if the kernel is too large, it escapes the domain, so the
	// normalization above must be corrected
	double m = 0;
	FORJ(h) FORI(w) m += g[j][i];
	FORJ(h) FORI(w) g[j][i] /= m;
}

static float average(float *x, int n)
{
	double r = 0;
	FORI(n) r += x[i]/n;
	return r;
}

// gaussian blur of a gray 2D image
static void gblur_gray(float *y, float *x, int w, int h, float s)
{
	s = 1/s;

	fftwf_complex *fx = fftwf_malloc(w*h*sizeof*fx);
	fft_2dfloat(fx, x, w, h);

	float *g = xmalloc(w*h*sizeof*g);
	fill_2d_gaussian_image(g, w, h, s);

	fftwf_complex *fg = fftwf_malloc(w*h*sizeof*fg);
	fft_2dfloat(fg, g, w, h);

	pointwise_complex_multiplication(fx, fx, fg, w*h);
	ifft_2dfloat(y, fx, w, h);

	fftwf_free(fx);
	fftwf_free(fg);
	free(g);
}

// gausian blur of a 2D image with pd-dimensional pixels
// (the blurring is performed independently for each co-ordinate)
static void gblur(float *y, float *x, int w, int h, int pd, float s)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h)
			c[i] = x[i*pd + l];
		gblur_gray(gc, c, w, h, s);
		FORI(w*h)
			y[i*pd + l] = gc[i];
	}
	free(c);
	free(gc);
}

#ifndef OMIT_GBLUR_MAIN
int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s s [in [out]]\n", *v);
		//                          0 1  2   3
		return EXIT_FAILURE;
	}
	float s = atof(v[1]);
	if (!s || !isfinite(s)) error("bad variance %g", s);
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	float *y = xmalloc(w*h*pd*sizeof*y);

	gblur(y, x, w, h, pd, s);

	iio_save_image_float_vec(out, y, w, h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
#endif//OMIT_GBLUR_MAIN
