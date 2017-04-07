#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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
	for (int i = 1; i < n; i++)
		if (x[i] > r)
			r = x[i];
	return r;
}

struct histogram {
	float binsize, binphase;
	int nbins; // always centered around zero
	int *count;
};

static void histogram_init(struct histogram *h,
		float binsize, float binphase, int nbins)
{
	h->binsize = binsize;
	h->binphase = binphase;
	h->nbins = nbins;
	h->count = malloc(nbins * sizeof*h->count);
	for (int i = 0; i < nbins; i++)
		h->count[i] = 0;
}

static void histogram_add(struct histogram *h, float x)
{
	int ix = h->nbins/2 + round( (x - h->binphase) / h->binsize );
	if (ix < 0) {
		ix = 0;
		fprintf(stderr, "histowarning (%g %g %d): %g out below\n",
				h->binsize, h->binphase, h->nbins, x);
	}
	if (ix >= h->nbins) {
		ix = h->nbins - 1;
		fprintf(stderr, "histowarning (%g %g %d): %g out above\n",
				h->binsize, h->binphase, h->nbins, x);
	}
	h->count[ix] += 1;
	//fprintf(stderr, "\tx=%g added to bin %d\n", x, ix);
}

// return the center of the bin with highest value
static float histogram_mode(struct histogram *h, bool solve_ties_upwards,
		int *count)
{
	int maxi;
	if (solve_ties_upwards) {
		maxi = h->nbins - 1;
		for (int i = 0; i < h->nbins; i++)
			if (h->count[i] > h->count[maxi])
				maxi = i;
	} else {
		maxi = 0;
		for (int i = h->nbins - 1; i >= 0; i--)
			if (h->count[i] > h->count[maxi])
				maxi = i;
	}
	if (count) *count = h->count[maxi];
	float fmaxi = (maxi - h->nbins/2) * h->binsize + h->binphase;
	return fmaxi;
}

static float histogram_debug(struct histogram *h)
{
	fprintf(stderr, "HISTOGRAM(bsize=%g bphase=%g nbins=%d)\n",
			h->binsize, h->binphase, h->nbins);
	for (int i = 0; i < h->nbins; i++)
		if (h->count[i])
			fprintf(stderr, "\th[%d] = %d\n", i, h->count[i]);
}

static float histomodev(float binsize, float binphase, float *x, int n,
		int *out_count)
{
	//fprintf(stderr, "binsize, binphase = %g %g\n", binsize, binphase);
	//fprintf(stderr, "we have %d numbers: ", n);
	//for (int i = 0; i < n; i++)
	//	fprintf(stderr, "%g%c", x[i], i==n-1?'\n':' ');
	float min = float_min(x, n);
	float max = float_max(x, n);
	//fprintf(stderr, "min, max = %g %g\n", min, max);
	int nbins=2+ceil(3*fmax(fabs(min+binphase),fabs(max+binphase))/binsize);
	//fprintf(stderr, "nbins = %d\n", nbins);
	struct histogram h[1];
	histogram_init(h, binsize, binphase, nbins);
	for (int i = 0; i < n; i++)
		histogram_add(h, x[i]);
	int rc;
	float r1 = histogram_mode(h, false, &rc);
	//float r2 = histogram_mode(h, true);
	//if (r1 != r2) fprintf(stderr, "WARNING: %g != %g (%g)\n",r1,r2,r2-r1);
	//histogram_debug(h);
	free(h->count);
	if(out_count) *out_count = rc;
	return r1;
}

#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *filename_out = pick_option(&c, &v, "i", "");
	if (c <= 4)
		return fprintf(stderr, "usage:\n\t"
				"%s binsize binphase x1 ... xN\n", *v);
	//                        0 1       2        3      N+3
	int n = c - 3;
	float x[n];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[3+i]);
	float binsize = atof(v[1]);
	float binphase = atof(v[2]);

	if (!*filename_out) {
		int rc;
		float r = histomodev(binsize, binphase, x, n, &rc);
		printf("%g (%d)\n", r, rc);
	} else {
		int w = 1024;
		float *t = malloc(2*w*w*sizeof*t);
		float first_phase = -binphase;
		float last_phase = binphase;
		float first_size = binsize/w;
		float last_size = binsize;
		fprintf(stderr, "phases(j) from %g to %g\n", first_phase, last_phase);
		fprintf(stderr, "sizes(i) from %g to %g\n", first_size, last_size);
		float first_freq = 1/first_size;
		float last_freq = 1/last_size;
		fprintf(stderr, "freq(i) from %g to %g\n", first_freq, last_freq);
		for (int j = 0; j < w; j++)
		for (int i = 0; i < w; i++)
		{
			float phase = first_phase+j*(last_phase-first_phase)/(w-1);
			float freq = first_freq+i*(last_freq-first_freq)/(w-1);
			int rc;
			t[2*(j*w+i)+0] = histomodev(1/freq, phase, x, n, &rc);
			t[2*(j*w+i)+1] = rc;
		}
		iio_write_image_float_vec(filename_out, t, w, w, 2);
	}

	return 0;
}
