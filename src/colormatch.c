#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// maximum pixel dimension for this program
#define MAX_DIM 100

// structure to describe the "contrast model" of an image
struct model_params {
	int pd;
	float mu[MAX_DIM];
	float sigma[MAX_DIM];
};

// structure to describe a contrast change between two images
struct model_map {
	int pd;
	float a[MAX_DIM];
	float b[MAX_DIM];
};

// apply a map to an image
static
void apply_map(float *y, struct model_map m, float *x, int w, int h, int pd)
{
	assert(m.pd == pd);
	for (int l = 0; l < pd; l++)
	for (int i = 0; i < w*h; i++)
		y[i*pd+l] = m.a[l] * x[i*pd+l] + m.b[l];
}

// function to compute a contrast change that maps one model to another
static
struct model_map match_models(struct model_params x, struct model_params y)
{
	assert(x.pd == y.pd);
	struct model_map m;
	m.pd = x.pd;
	for (int l = 0; l < m.pd; l++)
	{
		m.a[l] = x.sigma[l] / y.sigma[l];
		m.b[l] = x.mu[l] - m.a[l] * y.mu[l];
	}
	return m;
}

// compute avg and std of an array of numbers
static void compute_avg_std(float avg_std[2], float *x, int n)
{
	long double m = 0;
	for (int i = 0; i < n; i++)
		m += x[i];
	m /= n;

	long double s = 0;
	for (int i = 0; i < n; i++)
		s += (x[i] - m) * (x[i] - m);
	s = sqrt(s/n);

	avg_std[0] = m;
	avg_std[1] = s;
}

// computation of a model (very simple case)
static struct model_params compute_model(
		float *x, int w, int h, int pd,
		char *mask)
{
	struct model_params r;
	r.pd = pd;
	float *t = malloc(w*h*sizeof*t); // storage for valid samples
	for (int l = 0; l < pd; l++)
	{
		int cx = 0;
		for (int i = 0; i < w*h; i++)
			if (mask[i])
				t[cx++] = x[i*pd+l];
		float avg_std[2];
		compute_avg_std(avg_std, t, cx);
		r.mu[l]    = avg_std[0];
		r.sigma[l] = avg_std[1];
	}
	free(t);
	return r;
}


// extract a valid sample, nan otherwise
static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || j < 0 || i >= w || j >= h || l < 0 || l >= pd)
		return NAN;
	return x[(j*w+i)*pd+l];
}

#include "smapa.h"
SMART_PARAMETER_SILENT(VERBOSE,0)
static void dump_model_stderr(struct model_params p)
{
	for (int i = 0; i < p.pd; i++)
		fprintf(stderr, "p[%d] = %g %g\n", i, p.mu[i], p.sigma[i]);
}
static void dump_map_stderr(struct model_map m)
{
	for (int i = 0; i < m.pd; i++)
		fprintf(stderr, "m[%d] = %g %g\n", i, m.a[i], m.b[i]);
}

static void colormatch(float *C, float *A, int w, int h, int pd, float *B)
{
	// extract validity mask (pixel must be good on both images)
	char *mask = malloc(w*h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int good = 1;
		for (int l = 0; l < pd; l++)
		{
			if (!isfinite(getsample_nan(A, w, h, pd, i, j, l)))
				good = 0;
			if (!isfinite(getsample_nan(B, w, h, pd, i, j, l)))
				good = 0;
		}
		mask[j*w+i] = good;
	}

	// compute models of each image
	struct model_params p = compute_model(A, w, h, pd, mask);
	struct model_params q = compute_model(B, w, h, pd, mask);

	// compute mapping to transform one model into another
	struct model_map m = match_models(p, q);

	// apply transformation so that B looks like A
	apply_map(C, m, B, w, h, pd);

	if (VERBOSE() > 0)
	{
		dump_model_stderr(p);
		dump_model_stderr(q);
		dump_map_stderr(m);
	}
}

#include "iio.h"
int main(int c, char *v[])
{
	if (c != 4)
		return fprintf(stderr,"usage:\n\t%s in_A in_B out_BlikeA\n",*v);
		//                                0 1    2    3
	char *filename_in_A = v[1];
	char *filename_in_B = v[2];
	char *filename_out  = v[3];

	int w[2], h[2], pd[2];
	float *A = iio_read_image_float_vec(filename_in_A, w+0, h+0, pd+0);
	float *B = iio_read_image_float_vec(filename_in_B, w+1, h+1, pd+1);

	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "please use same size images\n");
	if (pd[0] != pd[1])
		return fprintf(stderr, "differnt pd match not implemented\n");

	float *C = malloc(w[1] * h[1] * *pd * sizeof*C);

	colormatch(C, A, *w, *h, *pd, B);

	iio_write_image_float_vec(filename_out, C, w[1], h[1], *pd);
	return 0;
}
