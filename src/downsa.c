#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

struct statistics_float {
	float min, max, median, average, sample, variance, middle, laverage;
};

#define STATISTIC_MEDIAN_BIAS 0
#define STATISTIC_MIDDLE_BIAS -1

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void statistics_getf_spoilable(struct statistics_float *s, float *f,
		int n)
{
	if (n == 0) {
		s->min = INFINITY;
		s->max = -INFINITY;
		s->median = NAN;
		s->average = 0;
		s->sample = NAN;
		s->variance = NAN;
		s->middle = NAN;
		s->laverage = NAN;
	}
	s->middle = f[n/2-1];
	int mt = STATISTIC_MEDIAN_BIAS;
	int mi = STATISTIC_MIDDLE_BIAS;
	switch(mi)
	{
		case -1: break;
		case 0: s->middle += f[n/2]; s->middle /=2; break;
		case 1: s->middle = f[n/2]; break;
		default: fail("bad STATISTIC_MEDIAN_BIAS %d", mt);
	}
	//
	qsort(f, n, sizeof*f, compare_floats);
	s->min = f[0];
	s->max = f[n-1];
	s->median = f[n/2-1];
	if (0 == n % 2)
	{
		int mtype = STATISTIC_MEDIAN_BIAS;
		switch(mtype)
		{
			case -1: break;
			case 0: s->median += f[n/2]; s->median /=2; break;
			case 1: s->median = f[n/2]; break;
			default: fail("bad STATISTIC_MEDIAN_BIAS %d", mtype);
		}
	}
	s->average = 0;
	for (int i = 0; i < n; i++)
		s->average += f[i];
	s->average /= n;
	s->laverage = 0;
	for (int i = 0; i < n; i++)
		s->laverage += exp(f[i]/255);
	s->laverage = log(s->laverage/n)*255;
	s->variance = 0;
	for (int i = 0; i < n; i++)
		s->variance = hypot(s->variance, s->average - f[i]);
	s->sample = f[randombounds(0, n-1)];
}

static void statistics_getf(struct statistics_float *s, float *fin, int n)
{
	if (n > 1000) {
		float *f = xmalloc(n*sizeof*f);
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
		xfree(f);
	} else {
		float f[n];
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
	}
}

static bool innerP(int w, int h, int i, int j)
{
	if (i < 0) return false;
	if (j < 0) return false;
	if (i >= w) return false;
	if (j >= h) return false;
	return true;
}

void downsa2d(float *oy, float *ox, int w, int h, int pd, int n, int ty)
{
	int W = w/n;
	int H = h/n;
	float (*x)[w][pd] = (void*)ox;
	float (*y)[W][pd] = (void*)oy;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	for (int l = 0; l < pd; l++)
	{
		int nv = 0;
		float vv[n*n];
		for (int jj = 0; jj < n; jj++)
		for (int ii = 0; ii < n; ii++)
		if (innerP(w, h, i*n+ii, j*n+jj))
		{
			float v = x[j*n+jj][i*n+ii][l];
			if (isfinite(v))
				vv[nv++] = v;
		}
		struct statistics_float s;
		statistics_getf(&s, vv, nv);
		float g;
		switch (ty)
		{
		case 'i': g = s.min;          break;
		case 'e': g = s.median;       break;
		case 'a': g = s.max;          break;
		case 'v': g = s.average;      break;
		case 'V': g = s.laverage;     break;
		case 's': g = s.variance;     break;
		case 'r': g = s.sample;       break;
		case 'n': g = nv;             break;
		case 'f': g = vv[0];          break;
		case 'l': g = vv[nv-1];       break;
		case 'c': g = vv[(nv-1)/2];   break;
		default:  fail("downsa type %c not implemented", ty);
		}
		y[j][i][l] = g;
	}
}

static char *help_string_name     = "downsa";
static char *help_string_version  = "downsa 1.0\n\nWritten by eml";
static char *help_string_oneliner = "zoom-out by combining NxN pixels into one";
static char *help_string_usage    = "usage:\n\t"
"downsa {i|e|a|v|r|n|f} N [in [out]]";
static char *help_string_long     =
"Downsa zooms-out an image by aggregation of NxN pixel values.\n"
"\n"
"Usage: downsa RULE n in out\n"
"   or: downsa RULE n in > out\n"
"   or: cat in | downsa RULE n > out\n"
"\n"
"Rules:\n"
" i     minimum\n"
" a     maximum\n"
" e     median\n"
" v     average\n"
" f     first\n"
" r     random choice among the NxN values\n"
" n     quantity of non-NANs\n"
"\n"
"Examples:\n"
" downsa a 2 stars.png small_stars.png              morphological zoom-out\n"
" cat big.png | blur g 1.6 | downsa f 2 > small     well-sampled zoom-out\n"
"\n"
"Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>.";
#include "help_stuff.c" // functions that print the strings named above

int main_downsa(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s {i|e|a|v|r|f} n [in [out]]\n", *v);
		//                         0        1      2  3   4
		return EXIT_FAILURE;
	}
	int n = atoi(v[2]);
	char *in = c > 3 ? v[3] : "-";
	char *out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	int W = w/n;
	int H = h/n;
	float *y = xmalloc(W*H*pd*sizeof*y);
	downsa2d(y, x, w, h, pd, n, v[1][0]);
	iio_write_image_float_vec(out, y, W, H, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_downsa(c, v); }
#endif
