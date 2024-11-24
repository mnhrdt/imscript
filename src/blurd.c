// blur of a 2D image, using various kernels

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>



#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
#endif

#include "fail.c"
#include "xmalloc.c"

static void *fftw_xmalloc(size_t n)
{
	double *r = fftw_malloc(n);
	if (!r)
		fail("coult not fftw_malloc %zu bytes\n", n);
	return r;
}

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORK(n) for(int k=0;k<(n);k++)
#define FORL(n) for(int l=0;l<(n);l++)



#ifndef USE_WISDOM
static void evoke_wisdom(void) {}
static void bequeath_wisdom(void) {}
#else//USE_WISDOM
#include "fftwisdom.c"
#endif//USE_WISDOM



// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2ddouble(fftw_complex *fx, double *x, int w, int h)
{
	fftw_complex *a = fftw_xmalloc(w*h*sizeof*a);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftw_plan p = fftw_plan_dft_2d(h, w, a, fx,
						FFTW_FORWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = x[i]; // complex assignment!
	fftw_execute(p);

	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_cleanup();
}

#include "smapa.h"
SMART_PARAMETER_SILENT(FIWARN,0)

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2ddouble(double *ifx,  fftw_complex *fx, int w, int h)
{
	fftw_complex *a = fftw_xmalloc(w*h*sizeof*a);
	fftw_complex *b = fftw_xmalloc(w*h*sizeof*b);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftw_plan p = fftw_plan_dft_2d(h, w, a, b,
						FFTW_BACKWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = fx[i];
	fftw_execute(p);
	double scale = 1.0/(w*h);
	FORI(w*h) {
		fftw_complex z = b[i] * scale;
		ifx[i] = crealf(z);
		if (FIWARN() > 0)
		{
			if (cimagf(z) > 0.001)
				fail("z is not real {cimagf(z)=%g} (set FIWARN=0 to run anyway)", cimagf(z));
			//assert(cimagf(z) < 0.001);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(a);
	fftw_free(b);
	fftw_cleanup();
}

SMART_PARAMETER_SILENT(BLUR_INVERSE,0)
SMART_PARAMETER_SILENT(BLUR_INVERSE_WIENER,0)
#define UGLY_HACK_FOR_WIENER_FILTERING 1

static void pointwise_complex_multiplication(fftw_complex *w,
		fftw_complex *z, fftw_complex *x, int n)
{
#ifdef UGLY_HACK_FOR_WIENER_FILTERING
	if (BLUR_INVERSE() > 0 || BLUR_INVERSE_WIENER() > 0) {
		if (BLUR_INVERSE_WIENER() > 0) {
			double t = BLUR_INVERSE_WIENER();
			FORI(n)
				w[i] = z[i] * x[i]/(cabs(x[i])*cabs(x[i])+t);
		} else
			FORI(n)
				w[i] = z[i] / x[i];
	}
	else
#endif//UGLY_HACK_FOR_WIENER_FILTERING
		FORI(n)
			w[i] = z[i] * x[i];
}

static void fill_2d_gaussian_image(double *gg, int w, int h, double inv_s)
{
	double (*g)[w] = (void *)gg;
	double alpha = inv_s*inv_s/(M_PI);

	FORJ(h) FORI(w) {
		double x = i < w/2 ? i : i - w;
		double y = j < h/2 ? j : j - h;
		double r = hypot(x, y);
		g[j][i] = alpha * exp(-r*r*inv_s*inv_s);
	}

	// if the kernel is too large, it escapes the domain, so the
	// normalization above must be corrected
	double m = 0;
	FORJ(h) FORI(w) m += g[j][i];
	FORJ(h) FORI(w) g[j][i] /= m;
}


// gaussian blur of a gray 2D image
void gblur_gray(double *y, double *x, int w, int h, double s)
{
	s = 1/s;

	fftw_complex *fx = fftw_xmalloc(w*h*sizeof*fx);
	fft_2ddouble(fx, x, w, h);

	double *g = xmalloc(w*h*sizeof*g);
	fill_2d_gaussian_image(g, w, h, s);

	fftw_complex *fg = fftw_xmalloc(w*h*sizeof*fg);
	fft_2ddouble(fg, g, w, h);

	pointwise_complex_multiplication(fx, fx, fg, w*h);
	ifft_2ddouble(y, fx, w, h);

	fftw_free(fx);
	fftw_free(fg);
	free(g);
}


// gausian blur of a 2D image with pd-dimensional pixels
// (the blurring is performed independently for each co-ordinate)
void gblur(double *y, double *x, int w, int h, int pd, double s)
{
	double *c = xmalloc(w*h*sizeof*c);
	double *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h) {
			double tmp = x[i*pd + l];
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


static double kernel_2d_square(double x, double y, double *p)
{
	int nx=0, ny=0;
	if (p[0] == 1)
		nx = ny = p[1];
	else if (p[0] == 2) {
		nx = p[1];
		ny = p[2];
	} else
		fail("square kernel needs 1 or 2 parameters (got %g)",p[0]);

	double r = 0;
	if (2*fabs(x) < nx && 2*fabs(y) < ny)
	       r = 1;

	return r;
}

static double kernel_2d_disk(double x, double y, double *p)
{
	double radius = p[1];

	double r = 0;
	if (hypot(x, y) < radius)
	       r = 1;

	return r;
}

static double kernel_2d_gaussian(double x, double y, double *p)
{
	double sigma = p[1];

	double a = x*x + y*y;
	double r = exp(-a/(2*sigma*sigma));
	return r;
}

static double kernel_2d_cauchy(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x,y)/sigma;
	double r = 1/(1+a*a);
	return r;
}

static double kernel_2d_bicauchy(double x, double y, double *p)
{
	double sigma = p[1];

	double a = (x*x + y*y)/(sigma*sigma);
	double r = 1/(1+a*a);
	return r;
}

static double kernel_2d_goodcauchy(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x,y)/sigma;
	double r = 1/pow(1+a*a, 1.5);
	return r;
}


static double kernel_2d_logcauchy(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x,y)/sigma;
	double r = a ? 1/(1+log(a)*log(a)) : 1;
	return r;
}

static double kernel_2d_powerlaw2(double x, double y, double *p)
{
	double sigma = p[1];

	double a = (x*x + y*y)/(sigma*sigma);
	double r = 1.0/(1.0 + a*a);
	return r;
}

static double kernel_2d_pareto(double x, double y, double *p)
{
	double alpha = p[1];

	double v = hypot(x, y);
	double r = v ? pow(v, alpha) : 1;
	return r;
}

static double kernel_2d_invr(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x, y) / sigma;
	double r = a ? 1/a : 1;
	return r;
}

static double kernel_2d_land(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x, y) / sigma;
	double r = a ? 1/(a*a*a) : 1;
	return r;
}

static double kernel_2d_riesz(double x, double y, double *p)
{
	double sigma = p[1];

	float n = hypot(x, y);
	float a = pow(hypot(x, y), sigma);
	float r = n ? pow(n, sigma-2) : 1;
	return r;
}

static double kernel_2d_ynvr(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x, y);
	double r = a ? 1/a : 1/sigma;
	return r;
}

SMART_PARAMETER_SILENT(LOGDESP,1.1)

static double kernel_2d_ilogr(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x, y) / sigma;
	double r = a ? 1/log(LOGDESP()+a) : 1;
	return r;
}

static double kernel_2d_logr(double x, double y, double *p)
{
	double sigma = p[1];
	double r = hypot(x, y);
	double v = r ? -log(r) : sigma;
	return v;
}


static double kernel_2d_r2logr(double x, double y, double *p)
{
	double sigma = p[1];

	double a = hypot(x, y) / sigma;
	double r = a ? a*a*log(a) : 1;
	return r;
}

static double kernel_2d_laplace(double x, double y, double *p)
{
	double sigma = p[1];

	double r = exp(-M_SQRT2*hypot(x,y)/sigma);
	return r;
}

static void fill_kernel_image(double *kk, int w, int h,
		double (*f)(double,double,double*), double *p)
{

	double (*k)[w] = (void *)kk;

	FORJ(h) FORI(w) {
		double x = i < w/2 ? i : i - w;
		double y = j < h/2 ? j : j - h;
		k[j][i] = f(x, y, p);
	}

	// if the kernel is too large, it escapes the domain, so the
	// normalization above must be corrected
	double m = 0;
	FORJ(h) FORI(w) m += k[j][i];
	FORJ(h) FORI(w) k[j][i] /= m;

}

static void substract_from_identity(double *k, int w, int h)
{
	for (int i = 0; i < w*h; i++)
		k[i] *= -1;
	k[0] += 1;
}

static void gray_fconvolution_2d(double *y, double *x, fftw_complex *fk,
		int w, int h)
{
	fftw_complex *fx = fftw_xmalloc(w*h*sizeof*fx);
	fft_2ddouble(fx, x, w, h);

	pointwise_complex_multiplication(fx, fx, fk, w*h);
	ifft_2ddouble(y, fx, w, h);

	fftw_free(fx);
}

static void color_fconvolution_2d(double *y, double *x, fftw_complex *fk,
		int w, int h, int pd)
{
	double *c = xmalloc(w*h*sizeof*c);
	double *kc = xmalloc(w*h*sizeof*kc);
	FORL(pd) {
		FORI(w*h) {
			double tmp = x[i*pd + l];
			if (!isfinite(tmp))
				tmp = 0;
			c[i] = tmp;//x[i*pd + l];
		}
		gray_fconvolution_2d(kc, c, fk, w, h);
		FORI(w*h)
			y[i*pd + l] = kc[i];
	}
	free(c);
	free(kc);
}

void blur_2d(double *y, double *x, int w, int h, int pd,
		char *kernel_id, double *param, int nparams)
{
	double p[1+nparams];
	for (int i = 0; i < nparams; i++)
		p[1+i] = param[i];
	p[0] = nparams;
	//FORI(nparams+1) fprintf(stderr, "p[%d:%d] = %g\n", i,nparams,p[i]);

	if (nparams == 1 && param[0] == 0) {
		for (int i = 0; i < w*h*pd; i++)
			y[i] = x[i];
		return;
	}

	double (*f)(double,double,double*) = NULL;
	switch(tolower(kernel_id[0])) {
	case 'g': f = kernel_2d_gaussian; break;
	case 'l': f = kernel_2d_laplace;  break;
	case 'c': f = kernel_2d_cauchy;   break;
	case 'k': f = kernel_2d_bicauchy;   break;
	case 'q': f = kernel_2d_logcauchy;   break;
	case 'u': f = kernel_2d_goodcauchy;   break;
	case 'd': f = kernel_2d_disk;     break;
	case 's': f = kernel_2d_square;   break;
	case 'p': f = kernel_2d_powerlaw2;   break;
	case 'a': f = kernel_2d_pareto;   break;
	case 'i': f = kernel_2d_invr;   break;
	case 'r': f = kernel_2d_riesz;   break;
	case 'y': f = kernel_2d_ynvr;   break;
	case 'z': f = kernel_2d_ilogr;   break;
	case 't': f = kernel_2d_r2logr;   break;
	case 'o': f = kernel_2d_logr;   break;
	default: fail("unrecognized kernel name \"%s\"", kernel_id);
	}

	double *k = xmalloc(w*h*sizeof*k);
	fill_kernel_image(k, w, h, f, p);
	if (isupper(kernel_id[0]))
		substract_from_identity(k, w, h);
	//void iio_write_image_double(char*,double*,int,int);
	//iio_write_image_double("/tmp/blurk.tiff", k, w, h);

	fftw_complex *fk = fftw_xmalloc(w*h*sizeof*fk);
	fft_2ddouble(fk, k, w, h);
	free(k);

	color_fconvolution_2d(y, x, fk, w, h, pd);

	fftw_free(fk);
}


#ifndef OMIT_BLUR_MAIN
#define MAIN_BLUR
#endif

#ifdef MAIN_BLUR

static char *help_string_name     = "blur";
static char *help_string_version  = "blur 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "smooth a 2D image using the selected kernel";
static char *help_string_usage    = "usage:\n\t"
"blur [-p 5] [-i] [-f]  [in [out]]";
static char *help_string_long     =
"Blur convolves the input image by the requested positive kernel.\n"
"Only the first letter of the kernel name is considered.\n"
"If the name of the kernel is uppercase, it subtracts the result\n"
"from the original image.\n"
"\n"
"Usage: blur KERNEL SIZE in.tiff out.tiff\n"
"   or: blur KERNEL SIZE in.tiff > out.tiff\n"
"   or: cat in.tiff | blur KERNEL SIZE > out.tiff\n"
"\n"
"Kernels:\n"
" square   a square block of the given radius\n"
" disk     a rasterized disk of the given radius\n"
" gauss    a Gaussian kernel of the given variance\n"
" laplace  a Laplace kernel of the given variance\n"
" cauchy   a Cauchy kernel of the given scale\n"
" q        Log-cauchy kernel\n"
" u        \"good-cauchy\"\n"
" p        powerlaw\n"
" a        pareto\n"
" i        inverse distance (useful for Shepard interpolation)\n"
" y        inverse distance (with different parameter normalization)\n"
" r        Riesz\n"
" z        inverse log-distance\n"
" t        r^2 log(r)  (useful for biharmonic interpolation)\n"
" o        log(r)\n"
"\n"
"Options:\n"
" -z        zero boundary\n"
" -s        symmetrized boundary\n"
" -p        periodic boundary\n"
"\n"
"Examples:\n"
" blur g 1.6                              Smooth an image by a slight amount\n"
" blur C 1 | qauto                        Linear retinex\n"
" plambda - \"x,l -1 *\" | blur i 0.25    Laplacian square root\n"
" plambda - \"x,l\" | blur z 0.25 | plambda - \"0 >\"      Linear dithering\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>.\n"
;
#include "help_stuff.c"
#include "parsenumbers.c"
#include "pickopt.c"
#include "iio.h"
int main_blur(int c, char *v[])
{
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	bool boundary_symmetric = pick_option(&c, &v, "s", NULL);
	bool boundary_periodic  = pick_option(&c, &v, "p", NULL);
	bool boundary_zero      = pick_option(&c, &v, "z", NULL);
	if (c != 5 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s kernel \"params\" [in [out]]\n", *v);
		//                0 1        2         3   4
		return EXIT_FAILURE;
	}
	char *kernel_id = v[1];
	char *kernel_params = v[2];
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int maxparam = 10;
	double param[maxparam];
	int nparams = parse_doubles(param, maxparam, kernel_params);
	if (nparams < 1) fail("please, give at least one parameter");
	//fprintf(stderr, "nparams = %d\n", nparams);
	//FORI(nparams)
	//	fprintf(stderr, "param[%d] = %g\n", i, param[i]);

	int w, h, pd;
	double *x = iio_read_image_double_vec(filename_in, &w, &h, &pd);
	double *y = xmalloc(w*h*pd*sizeof*y);

	if (boundary_symmetric || boundary_zero || boundary_periodic) {
		int ww = 2*w, hh = 2*h;
		double *xx = xmalloc(ww*hh*pd*sizeof*xx);
		double *yy = xmalloc(ww*hh*pd*sizeof*xx);
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int l = 0; l < pd; l++)
		{
			double g = x[(j*w+i)*pd+l];
			double g0 = boundary_zero ? 0 : g;
			int si = boundary_periodic ? w+i : 2*w - i - 1;
			int sj = boundary_periodic ? h+j : 2*h - j - 1;
			xx[( j*ww +  i)*pd+l] = g;
			xx[(sj*ww +  i)*pd+l] = g0;
			xx[( j*ww + si)*pd+l] = g0;
			xx[(sj*ww + si)*pd+l] = g0;
		}
		blur_2d(yy, xx, ww, hh, pd, kernel_id, param, nparams);
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int l = 0; l < pd; l++)
			y[(j*w+i)*pd+l] = yy[(j*ww+i)*pd+l];
		free(xx);
		free(yy);
	} else
		blur_2d(y, x, w, h, pd, kernel_id, param, nparams);

	iio_write_image_double_vec(filename_out, y, w, h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_blur(c, v); }
#endif

#endif//MAIN_BLUR
