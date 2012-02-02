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

#include "fail.c"
#include "xmalloc.c"


#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORK(n) for(int k=0;k<(n);k++)
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

// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued 3D image
static void fft_3dfloat(fftwf_complex *fx, float *x, int w, int h, int d)
{
	fftwf_complex *a = fftwf_malloc(w*h*d*sizeof*a);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_3d(d, h, w, a, fx,
						FFTW_FORWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h*d) a[i] = x[i]; // complex assignment!
	fftwf_execute(p);

	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_cleanup();
}

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial 3D image.
// The input data must be hermitic.
static void ifft_3dfloat(float *ifx,  fftwf_complex *fx, int w, int h, int d)
{
	fftwf_complex *a = fftwf_malloc(w*h*d*sizeof*a);
	fftwf_complex *b = fftwf_malloc(w*h*d*sizeof*b);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_3d(d, h, w, a, b,
						FFTW_BACKWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h*d) a[i] = fx[i];
	fftwf_execute(p);
	float scale = 1.0/(w*h*d);
	FORI(w*h*d) {
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

static void fill_3d_gaussian_image(float *gg, int w, int h, int d, float is[3])
{
	float (*g)[h][w] = (void *)gg;

	FORK(d) FORJ(h) FORI(w) {
		float x = i < w/2 ? i : i - w;
		float y = j < h/2 ? j : j - h;
		float z = k < d/2 ? k : k - d;
		float v = exp(-x*x*is[0] -y*y*is[1] -z*z*is[2]);
		g[k][j][i] = v;
	}

	double m = 0;
	FORK(d) FORJ(h) FORI(w) m += g[k][j][i];
	FORK(d) FORJ(h) FORI(w) g[k][j][i] /= m;
}

static float average(float *x, int n)
{
	double r = 0;
	FORI(n) r += x[i]/n;
	return r;
}

// gaussian blur of a gray 2D image
void gblur_gray(float *y, float *x, int w, int h, float s)
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

// gaussian blur of a gray 3D image
void gblur_gray_3d(float *y, float *x, int w, int h, int d, float rs[3])
{
	float s[3] = {1/rs[0], 1/rs[1], 1/rs[2]};
	int n = w * h * d;

	fftwf_complex *fx = fftwf_malloc(n*sizeof*fx);
	fft_3dfloat(fx, x, w, h, d);

	float *g = xmalloc(n*sizeof*g);
	fill_3d_gaussian_image(g, w, h, d, s);

	fftwf_complex *fg = fftwf_malloc(n*sizeof*fg);
	fft_3dfloat(fg, g, w, h, d);

	pointwise_complex_multiplication(fx, fx, fg, n);
	ifft_3dfloat(y, fx, w, h, d);

	fftwf_free(fx);
	fftwf_free(fg);
	free(g);
}

// gausian blur of a 2D image with pd-dimensional pixels
// (the blurring is performed independently for each co-ordinate)
void gblur(float *y, float *x, int w, int h, int pd, float s)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h) {
			float tmp = x[i*pd + l];
			if (!isfinite(tmp))
				tmp = 0;
			c[i] = tmp;//x[i*pd + l];
		}
		if (s)
			gblur_gray(gc, c, w, h, s);
		else FORI(w*h) gc[i] = c[i];
		FORI(w*h)
			y[i*pd + l] = gc[i];
	}
	free(c);
	free(gc);
}

// gausian blur of a 3D image with pd-dimensional pixels
// (the blurring is performed independently for each co-ordinate)
void gblur3d(float *y, float *x, int w, int h, int d, int pd, float s[3])
{
	int n = w * h * d;
	float *c = xmalloc(n*sizeof*c);
	float *gc = xmalloc(n*sizeof*gc);
#pragma omp parallel for
	FORL(pd) {
		FORI(n) {
			float tmp = x[i*pd + l];
			if (!isfinite(tmp))
				tmp = 0;
			c[i] = tmp;//x[i*pd + l];
		}
		gblur_gray_3d(gc, c, w, h, d, s);
		FORI(n)
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
	if (!s || !isfinite(s)) fail("bad variance %g", s);
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
