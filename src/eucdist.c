// euclidean distance transform (Mejister-Felzenszwalb-Huttenlocher algorithm)

// this implementation is based on the article
// "Distance Transforms of Sampled Functions",
// by P.F.Felzenszwalb and D.P.Huttenlocher,
// published on "Theory of Computing" in 2012

#include <math.h>

static float obtain_slope(float *f, int q, int p)
{
	if (!isfinite(f[q]) && !isfinite(f[p])) return 0;
	return ((f[q] + q*q) - (f[p] + p*p)) / (2*q - 2*p);
}

// algorithm 1 from the paper, using the same variable names
static void squared_distances_1d(
		float *d,  // output distance
		float *f,  // input data
		int n      // length
		)
{
	float z[n+1];  // temporary array
	int v[n];      // index table
	int k = 0;     // last index

	v[0] = 0;
	z[0] = -INFINITY;
	z[1] =  INFINITY;
	for (int q = 1; q < n; q++) {
		// TODO: simplify the following ugly loop
		float s;
		while ((s = obtain_slope(f, q, v[k])) < z[k])
			k--;
		k++;
		v[k] = q;
		z[k] = s;
		z[k+1] = INFINITY;
	}

	k = 0;
	for (int q = 0; q < n; q++) {
		while (z[k+1] < q)
			k++;
		d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
	}
}

static void squared_distances_2d(float *f, int w, int h)
{
	// columns
	for (int i = 0; i < w; i++)
	{
		float t[h], d[h];
		for (int j = 0; j < h; j++) t[j] = f[j*w+i];
		squared_distances_1d(d, t, h);
		for (int j = 0; j < h; j++) f[j*w+i] = d[j];
	}

	// rows
	for (int j = 0; j < h; j++)
	{
		float t[w], d[w];
		for (int i = 0; i < w; i++) t[i] = f[j*w+i];
		squared_distances_1d(d, t, w);
		for (int i = 0; i < w; i++) f[j*w+i] = d[i];
	}
}

void squared_euclidean_distance_to_nonzeros(
		float *x,   // input/output: (binary image/distance to nonzeros)
		int w,      // width
		int h       // height
		)
{
	for (int i = 0; i < w*h; i++)
		x[i] = x[i] > 0 ? 0: INFINITY;
	squared_distances_2d(x, w, h);
}

#include <stdlib.h>
void signed_distance_to_mask(float *x, int w, int h)
{
	float *tp = malloc(w * h * sizeof*tp);
	float *tm = malloc(w * h * sizeof*tm);
	for (long i = 0; i < w*h; i++)
	{
		tp[i] = x[i] > 0;
		tm[i] = !tp[i];
	}
	squared_euclidean_distance_to_nonzeros(tp, w, h);
	squared_euclidean_distance_to_nonzeros(tm, w, h);
	for (long i = 0; i < w*h; i++)
		x[i] = x[i] > 0 ? sqrt(tm[i])-0.5:0.5-sqrt(tp[i]);
	free(tp);
	free(tm);
}

#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	_Bool s = pick_option(&c, &v, "s", NULL);
	if (c != 3) return 3;
	int w, h;
	float *x = iio_read_image_float(v[1], &w, &h);
	if (!s)
		squared_euclidean_distance_to_nonzeros(x, w, h);
	else
		signed_distance_to_mask(x, w, h);

	iio_write_image_float(v[2], x, w, h);
	return 0;
}
