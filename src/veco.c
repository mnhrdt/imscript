#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <math.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

static float float_sum(float *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return r;
}

static float float_avg(float *x, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i];
	return n?r/n:r;
}

static float float_mul(float *x, int n)
{
	double r = 1;
	for (int i = 0; i < n; i++)
		r *= x[i];
	return r;
}

static float float_min(float *x, int n)
{
	float r = INFINITY;
	for (int i = 0; i < n; i++)
		if (x[i] < r)
			r = x[i];
	return r;
}

static float float_max(float *x, int n)
{
	float r = -INFINITY;
	for (int i = 0; i < n; i++)
		if (x[i] > r)
			r = x[i];
	return r;
}

int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static float float_med(float *x, int n)
{
	if (!n) return NAN;//fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return x[0];
	qsort(x, n, sizeof*x, compare_floats);
	return x[n/2];
}

static float float_cnt(float *x, int n)
{
	return n;
}

static float float_logavg(float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += log(x[i]);
	return n?exp(r/n):0;
}

static float float_logsumexp(float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += exp(x[i]);
	return n?log(r/n):0;
}

static float float_std(float *x, int n)
{
	float m = float_avg(x, n);
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += (x[i] - m) * (x[i] - m);
	return sqrt(r/n);
}

static float float_iqd(float *x, int n)
{
	if (!n) return NAN;//fail("empty list of pixel values!");
	if (n == 1) return 0;
	qsort(x, n, sizeof*x, compare_floats);
	if (n == 2) return x[1] - x[0];
	if (n == 3) return x[2] - x[0];
	if (n == 4) return x[2] - x[1];
	int a = round(0.25 * (n-1));
	int b = round(0.75 * (n-1));
	return x[b] - x[a];
}

static float float_mod(float *x, int n)
{
	float h[0x100];
	for (int i = 0; i < 0x100; i++)
		h[i] = 0;
	for (int i = 0; i < n; i++)
	{
		int xi = x[i];
		if (xi < 0) continue;//fail("negative xi=%g", x[i]);//xi = 0;
		if (xi > 0xff) continue;//fail("large xi=%g", x[i]);//xi = 0xff;
		h[xi] += 2;
		if (xi > 0) h[xi-1] += 1;
		if (xi < 0xff) h[xi+1] += 1;
	}
	int mi = 0x80;
	for (int i = 0; i < 0x100; i++)
		if (h[i] > h[mi])
			mi = i;
	return mi;
}

static float float_modH(float *x, int n)
{
	static float p = 1;
	if (n == -1)
		return p = *x;
	float h[0x100];
	for (int i = 0; i < 0x100; i++)
		h[i] = 0;
	for (int i = 0; i < n; i++)
	{
		int xi = p*floor(x[i]/p);
		if (xi < 0) continue;//fail("negative xi=%g", x[i]);//xi = 0;
		if (xi > 0xff) continue;//fail("large xi=%g", x[i]);//xi = 0xff;
		h[xi] += 2;
		if (xi > 0) h[xi-1] += 1;
		if (xi < 0xff) h[xi+1] += 1;
	}
	int mi = 0x80;
	for (int i = 0; i < 0x100; i++)
		if (h[i] > h[mi])
			mi = i;
	return mi;
}

static float float_harmonic(float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += 1/x[i];
	return n/r;
}

static float float_euclidean(float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r = hypot(r, x[i]);
	return r/sqrt(n);
}

static float float_geometric(float *x, int n)
{
	long double r = 1;
	for (int i = 0; i < n; i++)
		r *= x[i];
	return powl(r, 1.0/n);
}

static float float_holder(float *x, int n)
{
	static float p = 1;
	if (n == -1)
		return p = *x;
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += powl(x[i], p);
	return powl(r/n, 1/p);
}

static float float_lehmer(float *x, int n)
{
	static float p = 2;
	if (n == -1)
		return p = *x;
	long double a = 0;
	long double b = 0;
	for (int i = 0; i < n; i++)
	{
		a += powl(x[i], p);
		b += powl(x[i], p-1);
	}
	return a / b;
}

static float float_gini(float *x, int n)
{
	static float p = 3;
	static float q = 2;
	if (n == -1) {
		p = x[0];
		q = x[1];
		return p*q;
	}
	if (p == q) {
		//long double a = 1;
		//long double b = 0;
		//for (int i = 0; i < n; i++)
		//{
		//	long double xip = powl(x[i], p);
		//	b += xip;
		//	a *= powl(x[i], xip);
		//}
		//return powl(a, 1/b);
		//
		// the following code is algebraically equivalent,
		// but much better conditionated than the formula above
		long double a = 0;
		long double b = 0;
		for (int i = 0; i < n; i++)
		{
			long double xip = powl(x[i], p);
			a += xip * logl(x[i]);
			b += xip;
		}
		return exp(a/b);
	} else {
		long double a = 0;
		long double b = 0;
		for (int i = 0; i < n; i++)
		{
			a += powl(x[i], p);
			b += powl(x[i], q);
		}
		return powl(a / b, 1 / (p - q));
	}
}

#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

static float float_medv(float *x, int n)
{
	if (!n) fail("empty list of pixel values!");
	if (n == 1) return x[0];
	if (n == 2) return (x[0] + x[1])/2;
	qsort(x, n, sizeof*x, compare_floats);
	if (EVENP(n))
		return (x[n/2] + x[-1+n/2])/2;
	else
		return x[n/2];
}

static float float_percentile(float *x, int n)
{
	static float percentile = 50;
	if (n == -1) { percentile = fmax(0, fmin(*x, 100)); return 0; }
	if (n == 0) return NAN;
	if (n == 1) return x[0];
	qsort(x, n, sizeof*x, compare_floats);
	int i = round(percentile * (n - 1) / 100.0);
	//fprintf(stderr, "n=%d, percentile=%g, i=%d\n", n, percentile, i);
	assert(i >= 0);
	assert(i < n);
	return x[i];
}

static float EEE(float *x, int n, float p, float m)
{
	long double r = 0;
	for (int i = 0; i < n; i++)
		r += powl(fabs(x[i] - m),  p);
	return r;
}

static float float_pargmineg(float *x, float p, int n)
{
	long double score[n];

	for (int i = 0; i < n; i++)
		score[i] = 0;

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	if (j != i)
		score[i] += powl(fabs(x[i] - x[j]), p);

	int ridx = 0;
	if (p > 0) {
		for (int i = 1; i < n; i++)
		if (score[i] < score[ridx])
			ridx = i;
	}
	if (p < 0) {
		for (int i = 1; i < n; i++)
		if (score[i] > score[ridx])
			ridx = i;
	}

	return x[ridx];
}

static float float_pargmin(float *x, int n)
{
	static float p = 2;
	if (n == -1)
		return p = *x;
	if (p < 1) return float_pargmineg(x, p, n);
	if (p == 0) return NAN;
	long double mo = 0;
	long double m = float_avg(x, n);
	for (int j = 0; j < 25; j++)
	{
		fprintf(stderr, "m = %g     E(m,%g)=%g\n", (double)m, p, EEE(x,n,p,m));
		long double a = 0;
		long double b = 0;
		for (int i = 0; i < n; i++)
		{
			long double w = powl(fabs(x[i] - m), p - 2);
			fprintf(stderr, "\tw = |%g - %g|^%g = %g\n", x[i], (double)m, p-2, (double)w);
			a += w * x[i];
			b += w;
		}
		m = a / b;
		if (fabs(mo - m) < fabs(m)*1e-6) break;
		mo = m;
	}
	if (1) {
		float mt[5];
		mt[0] = float_min(x, n);
		mt[1] = float_avg(x, n);
		mt[2] = float_max(x, n);
		mt[3] = float_medv(x, n);
		mt[4] = (mt[0] + mt[2])/2;
		fprintf(stderr,"E_%g(min=\t%g)\t%g\n",p,mt[0],EEE(x,n,p,mt[0]));
		fprintf(stderr,"E_%g(avg=\t%g)\t%g\n",p,mt[1],EEE(x,n,p,mt[1]));
		fprintf(stderr,"E_%g(max=\t%g)\t%g\n",p,mt[2],EEE(x,n,p,mt[2]));
		fprintf(stderr,"E_%g(med=\t%g)\t%g\n",p,mt[3],EEE(x,n,p,mt[3]));
		fprintf(stderr,"E_%g(cen=\t%g)\t%g\n",p,mt[4],EEE(x,n,p,mt[4]));
		fprintf(stderr,"E_%g( m =\t%g)\t%g\n",p,(double)m,EEE(x,n,p,m));
	}
	fprintf(stderr, " = %g     E(m,%g)=%g\n", (double)m, p, EEE(x,n,p,m));
	return m;
}


static float float_first(float *x, int n)
{
	if (n)
		return x[0];
	else
		fail("empty list of pixel values!");
}

static float float_pick(float *x, int n)
{
	if (n) {
		int i = randombounds(0, n-1);
		return x[i];
	}
	else
		fail("empty list of pixel values!");
}

typedef bool (*isgood_t)(float);

static bool isgood_finite(float x) { return isfinite(x); }
static bool isgood_always(float x) { return true; }
static bool isgood_numeric(float x) { return !isnan(x); }
static bool isgood_nonzero(float x) { return x != 0; }
static bool isgood_positive(float x) { return x >= 0; }
static bool isgood_negative(float x) { return x < 0; }

static char *help_string_name     = "veco";
static char *help_string_version  = "veco 1.0\n\nWritten by eml";
static char *help_string_oneliner = "combine several scalar images into one";
static char *help_string_usage    = "usage:\n\t"
"veco {sum|min|max|mul|med|...} in1 in2 ... {> out|-o out}";
static char *help_string_long     =
"Veco combines several scalar images by a pixelwise operation\n"
"\n"
"Usage: veco OPERATION in1 in2 in3 ... > out\n"
"   or: veco OPERATION in1 in2 in3 ... -o out\n"
"\n"
"Options:\n"
" -o file      use a named output file instead of stdout\n"
" -g GOODNESS  use a specific GOODNESS criterion for discarding samples\n"
" -x NUMBERS   instead of images, combine the given numbers\n"
" -k IMAGE     operate over the channels of a single image\n"
" -c IMAGE     operate over the columns of a single image\n"
" -i           operate independently along the dimensions of multispectral images\n"
"\n"
"Operations:\n"
" min          minimum value of the good samples\n"
" max          maximum value of the good samples\n"
" avg          average of the good samples\n"
" sum          sum of the good samples\n"
" med          medoid of the good samples (the central sample, rounding down)\n"
" medv         median of the good samples (average of 1 or 2 central samples)\n"
" mod          mode of the good samples (the value that appears more times)\n"
" cnt          number of good samples\n"
" mul          product of all good samples\n"
" first        the first good sample\n"
" rnd          a randomly chosen good sample\n"
" qX           Xth percentile\n"
" MP           Pth power mean\n"
" FP           Fréchet typical position, or minimizer of P-norm (M1=V2, etc)\n"
" LP           P-Lehmer mean\n"
" GP,Q         Gini's mean with parameters P and Q\n"
" euc          euclidean norm (M2)\n"
" geo          geometric mean (M0)\n"
" har          harmonic mean (M-1)\n"
" lav          logarithmic average\n"
" lse          log-sum-exp (a.k.a. soft max)\n"
" std          standard deviation\n"
" iqd          interquartile distance\n"
"\n"
"Goodness criteria:\n"
" numeric      whether the sample is not NAN, this is the default\n"
" finite       whether the sample is a finite number\n"
" always       consider all samples regardless of their value\n"
" nonzero      whether the sample is numeric and not 0 or -0\n"
" positive     whether the sample is >= 0\n"
"\n"
"Examples:\n"
" veco avg i*.png -o avg.png     Compute the average of a bunch of images\n"
" veco M0 -x 1 2 3               Compute the geometric mean of three numbers\n"
"\n"
"Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>."
;
#include "help_stuff.c"
#include "pickopt.c"
int main_veco(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);
	bool use_numbers =   pick_option(&c, &v, "x", 0);
	bool by_columns  =   pick_option(&c, &v, "c", 0);
	bool by_channels =   pick_option(&c, &v, "k", 0);
	bool by_full_img =   pick_option(&c, &v, "f", 0);
	bool indep_chans =   pick_option(&c, &v, "i", 0);
	char *goodness   =   pick_option(&c, &v, "g", "numeric");
	char *filename_out = pick_option(&c, &v, "o", "-");
	if (c < 3 && !by_channels) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg|mul|med} [v1 ...] > out\n", *v);
		//          0  1                          2  3
		return EXIT_FAILURE;
	}
	int n = c - 2;
	char *operation_name = v[1];
	float (*f)(float *,int) = NULL;
	if (0 == strcmp(operation_name, "sum"))   f = float_sum;
	if (0 == strcmp(operation_name, "mul"))   f = float_mul;
	if (0 == strcmp(operation_name, "prod"))  f = float_mul;
	if (0 == strcmp(operation_name, "avg"))   f = float_avg;
	if (0 == strcmp(operation_name, "min"))   f = float_min;
	if (0 == strcmp(operation_name, "max"))   f = float_max;
	if (0 == strcmp(operation_name, "med"))   f = float_med;
	if (0 == strcmp(operation_name, "mod"))   f = float_mod;
	if (0 == strcmp(operation_name, "cnt"))   f = float_cnt;
	if (0 == strcmp(operation_name, "medi"))   f = float_med;
	if (0 == strcmp(operation_name, "medv"))   f = float_medv;
	if (0 == strcmp(operation_name, "rnd"))   f = float_pick;
	if (0 == strcmp(operation_name, "euc"))   f = float_euclidean;
	if (0 == strcmp(operation_name, "geo"))   f = float_geometric;
	if (0 == strcmp(operation_name, "har"))   f = float_harmonic;
	if (0 == strcmp(operation_name, "lav"))   f = float_logavg;
	if (0 == strcmp(operation_name, "lse"))   f = float_logsumexp;
	if (0 == strcmp(operation_name, "std"))   f = float_std;
	if (0 == strcmp(operation_name, "iqd"))   f = float_iqd;
	if (0 == strcmp(operation_name, "first")) f = float_first;
	if (*operation_name == 'q') {
		float p = atof(1 + operation_name);
		f = float_percentile;
		f(&p, -1);
	}
	if (*operation_name == 'M') {
		float p = atof(1 + operation_name);
		f = float_holder;
		f(&p, -1);
		if (p == 0)         f = float_geometric;
		if (p == INFINITY)  f = float_max;
		if (p == -INFINITY) f = float_min;
	}
	if (*operation_name == 'F') {
		float p = atof(1 + operation_name);
		f = float_pargmin;
		f(&p, -1);
	}
	if (*operation_name == 'H') {
		float p = atof(1 + operation_name);
		f = float_modH;
		f(&p, -1);
	}
	if (*operation_name == 'L') {
		float p = atof(1 + operation_name);
		f = float_lehmer;
		f(&p, -1);
	}
	if (*operation_name == 'G') {
		float pq[2];
		int r = sscanf(operation_name, "G%g,%g", pq, pq+1);
		if (r != 2) fail("unrecognized op \"%s\"", operation_name);
		f = float_gini;
		f(pq, -1);
	}
	if (!f) fail("unrecognized operation \"%s\"", operation_name);
	bool (*isgood)(float) = NULL;
	if (0 == strcmp(goodness, "finite"))   isgood = isgood_finite;
	if (0 == strcmp(goodness, "numeric"))  isgood = isgood_numeric;
	if (0 == strcmp(goodness, "always"))   isgood = isgood_always;
	if (0 == strcmp(goodness, "nonzero"))  isgood = isgood_nonzero;
	if (0 == strcmp(goodness, "positive")) isgood = isgood_positive;
	if (!isgood) fail("unrecognized goodness \"%s\"", goodness);
	if (use_numbers) {
		float x[n];
		for (int i = 0; i < n; i++)
			x[i] = atof(v[i+2]);
		float tmp[n];
		int ngood = 0;
		for (int i = 0; i < n; i++)
			if (isgood(x[i]))
				tmp[ngood++] = x[i];
		float y = f(tmp, ngood);
		printf("%g\n", y);
	} else if (n == 1 && by_columns) {
		int w, h;
		float *x = iio_read_image_float(v[2], &w, &h);
		float *y = xmalloc(w * sizeof*y);
		for (int i = 0; i < w; i++)
		{
			float tmp[h];
			int ngood = 0;
			for (int j = 0; j < h; j++)
				if (isgood(x[j*w+i]))
					tmp[ngood++] = x[j*w+i];
			y[i] = f(tmp, ngood);
		}
		iio_write_image_float(filename_out, y, w, 1);
	} else if (n < 2 && by_channels) {
		int w, h, pd;
		float *x = iio_read_image_float_vec(c>2?v[2]:"-", &w, &h, &pd);
		float *y = xmalloc(w * h * sizeof*y);
		for (int i = 0; i < w*h; i++)
		{
			float tmp[pd];
			int ngood = 0;
			for (int j = 0; j < pd; j++)
				if (isgood(x[i*pd+j]))
					tmp[ngood++] = x[i*pd+j];
			y[i] = f(tmp, ngood);
		}
		iio_write_image_float(filename_out, y, w, h);
	} else if (n < 2 && by_full_img) {
		int w, h, pd;
		float *x = iio_read_image_float_vec(c>2?v[2]:"-", &w, &h, &pd);
		float *y = xmalloc(w * h * sizeof*y);
		int ngood = 0;
		for (int i = 0; i < w*h; i++)
		{
			for (int j = 0; j < pd; j++)
				if (isgood(x[i*pd+j]))
					y[ngood++] = x[i*pd+j];
		}
		float out = f(y, ngood);
		printf("%lf\n", out);
	} else if (indep_chans) {
		float *x[n];
		int w[n], h[n], d[n];
		for (int i = 0; i < n; i++)
			x[i] = iio_read_image_float_split(v[i+2], w+i,h+i,d+i);
		for (int i = 0; i < n; i++) {
			if (w[i] != *w || h[i] != *h || d[i] != *d)
				fail("%dth image size mismatch\n", i);
		}
		float (*y) = xmalloc(*w * *h * *d * sizeof*y);
		for (int l = 0; l < *d; l++)
		{
			int O = l * *w * *h;
			for (int i = 0; i < *w * *h; i++)
			{
				float tmp[n];
				int ngood = 0;
				for (int j = 0; j < n; j++)
					if (isgood(x[j][i+O]))
						tmp[ngood++] = x[j][i+O];
				y[i+O] = f(tmp, ngood);
			}
		}
		iio_write_image_float_split(filename_out, y, *w, *h, *d);
	} else {
		float *x[n];
		int w[n], h[n];
		for (int i = 0; i < n; i++)
			x[i] = iio_read_image_float(v[i+2], w + i, h + i);
		for (int i = 0; i < n; i++) {
			if (w[i] != *w || h[i] != *h)
				fail("%dth image size mismatch\n", i);
		}
		float (*y) = xmalloc(*w * *h * sizeof*y);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < *w * *h; i++)
		{
			float tmp[n];
			int ngood = 0;
			for (int j = 0; j < n; j++)
				if (isgood(x[j][i]))
					tmp[ngood++] = x[j][i];
			y[i] = f(tmp, ngood);
		}
		iio_write_image_float(filename_out, y, *w, *h);
	}
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_veco(c, v); }
#endif
