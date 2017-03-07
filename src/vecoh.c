#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

#include "smapa.h"
SMART_PARAMETER(PRECISION,10)
SMART_PARAMETER_SILENT(NUMBER_OF_KMEANS_ITERATIONS,5)


int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

// compute the avg of the variances of each cluster
static float average_of_cluster_variances(
		int *group, float *mean, int k,
		float *x, int n)
{
	if (!k) return 0;
	float var[k], cnt[k];
	for (int i = 0; i < k; i++)
		var[i] = cnt[i] = 0;

	for (int i = 0; i < n; i++)
	{
		float vi = x[i];
		int gi = group[i];
		assert(gi >= 0);
		assert(gi <  k);
		var[gi] += (x[i] - mean[gi]) * (x[i] - mean[gi]);
		cnt[gi] += 1;
	}
	for (int i = 0; i < k; i++)
		var[i] /= cnt[i];

	float r = 0;
	for (int i = 0; i < k; i++)
		r += var[i];
	return r / k;
}

// given a set of centroids,
// assign the index of the nearest centroid to each point
static void assign_each_point_to_nearest_mean(int *assign, float *mean,
		int k, float *x, int n)
{
	for (int i = 0; i < n; i++)
	{
		int bj = 0;
		for (int j = 0; j < k; j++)
			if (fabs(x[i]-mean[j]) < fabs(x[i]-mean[bj]))
				bj = j;
		assign[i] = bj;
	}
}

// the standard k-means algorithm for 1d points
static void kmeans_1d(float *means, float *x, int n, int k)
{
	// order and sanity
	for (int i = 0; i < k; i++)
		means[i] = NAN;
	if (!n) return;
	qsort(x, n, sizeof*x, compare_floats);
	if (n <= k) {
		for (int i = 0; i < n; i++)
			means[i] = x[i];
		return;
	}

	// initialize the position of the clusters uniformly
	assert(n > k);
	for (int i = 0; i < k; i++)
	{
		int idx = lrint(  (i + 1.0) * (n - 1.0) / (k + 1.0)  );
		assert(idx >= 0);
		assert(idx < n);
		means[i] = x[idx];
	}
	for (int i = 1; i < k; i++)
		if (means[i-1] == means[i])
		{
			means[i-1] -= 0.01;
			means[i] += 0.01;
		}

	// run a few iterations of k-means
	int numit = NUMBER_OF_KMEANS_ITERATIONS();
	for (int cx = 0; cx < numit; cx++)
	{
		// assign each point to the nearest mean
		int assignement[n];
		for (int i = 0; i < n; i++)
		{
			int bj = 0;
			for (int j = 0; j < k; j++)
				if (fabs(x[i]-means[j]) < fabs(x[i]-means[bj]))
					bj = j;
			assignement[i] = bj;
		}

		// update each mean as the mean of the assigned points
		int count[k];
		for (int j = 0; j < k; j++)
			means[j] = count[j] = 0;
		for (int i = 0; i < n; i++)
		{
			means[assignement[i]] += x[i];
			count[assignement[i]] += 1;
		}
		for (int j = 0; j < k; j++)
			if (count[j])
				means[j] /= count[j];
			else
				means[j] = NAN;
	}

}

// INPUT: n samples x[i], desired output dimension "dy"
// OUTPUT:
// y[0] = number of modes of the samples x[i] (must satisfy y[0] < dy-1)
// y[1] = position of first mode
// y[2] = position of second mode
// ...
// y[y[0]] = position of the last mode
// y[y[0]+1] = nan
// ...
// y[dy-1] = nan
static void float_kmeans(float *y, int dim_y, float *x, int n, int k)
{
	for (int i = 0; i < dim_y; i++) y[i] = NAN;
	if (k >= dim_y) k = dim_y - 1;
	float means[k];
	kmeans_1d(means, x, n, k);
	int count = 0;
	for (int i = 0; i < k; i++)
		count += isfinite(means[i]);
	y[0] = count;
	int ypos = 1;
	for (int i = 0; i < k; i++)
		if (isfinite(means[i]))
			y[ypos++] = means[i];
}

// run k-means with increasing k until the average variance of the clusters
// falls below a pre-sed precision
static void float_varkmeans(float *y, int dim_y, float *x, int n, float prec)
{
	for (int pre_k = 1; pre_k < n; pre_k++)
	{
		float_kmeans(y, dim_y, x, n, pre_k);
		int k = y[0];
		int group[n];
		assign_each_point_to_nearest_mean(group, y+1, k, x, n);
		float acv = average_of_cluster_variances(group, y+1, k, x, n);
		if (acv < prec*prec)
			break;
	}
}


// main for debugging: try the k-means clustering for k=env(KKK)
SMART_PARAMETER(KKK,2)
static int main_ttry(int c, char *v[])
{
	if(c<3)return fprintf(stderr,"usage:\n\tkmeans x1 ... xn\n");
	int n = c - 1;
	int max_k = 20;
	float x[n], y[max_k];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[i+1]);
	int kkk = KKK();
	float_kmeans(y, max_k, x, n, kkk);
	int k = y[0];
	int group[n];
	assign_each_point_to_nearest_mean(group, y+1, k, x, n);
	float acv = average_of_cluster_variances(group, y+1, k, x, n);
	fprintf(stderr, "racv = %g\t", sqrt(acv));
	fprintf(stderr, "got %d centers (", k);
	for (int i = 0; i < k; i++)
		fprintf(stderr, "%g%c", y[1+i], i+1==k?')':' ');
	puts("");
	return 0;
}

// main for debugging: try the k-estimation a given PRECISION
static int main_xtry(int c, char *v[])
{
	if(c<3)return fprintf(stderr,"usage:\n\tkmeans x1 ... xn\n");
	int n = c - 1;
	int max_k = 20;
	float x[n], y[max_k];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[i+1]);
	float_varkmeans(y, max_k, x, n, PRECISION());
	int k = y[0];
	int group[n];
	assign_each_point_to_nearest_mean(group, y+1, k, x, n);
	float acv = average_of_cluster_variances(group, y+1, k, x, n);
	fprintf(stderr, "racv = %g\t", sqrt(acv));
	fprintf(stderr, "got %d centers (", k);
	for (int i = 0; i < k; i++)
		fprintf(stderr, "%g%c", y[1+i], i+1==k?')':' ');
	puts("");
	return 0;
}

// wether a floating point vector is good or not
static bool isgood(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isfinite(x[i]))
			return false;
	return true;
}

#include "pickopt.c"
SMART_PARAMETER_SILENT(VECOH_VERBOSE,0)
// main main: compute the modes of a series of images
int main_vecoh(int c, char *v[])
{
	if (pick_option(&c, &v, "t", 0)) return main_ttry(c, v);
	if (pick_option(&c, &v, "x", 0)) return main_xtry(c, v);
	char *filename_out = pick_option(&c, &v, "o", "-");
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {kmeans|modes} [v1 ...] [-o out]\n", *v);
		//          0  1              2  3
		return 1;
	}
	int n = c - 2;
	char *operation_name = v[1];
	//void (*f)(float*,float*,int,int) = NULL;
	//if (0 == strcmp(operation_name, "kmeans")) f = float_kmeans;
	//if (0 == strcmp(operation_name, "modes"))  f = float_modes;
	//if (!f) fail("unrecognized operation \"%s\"", operation_name);
	float *x[n];
	int w[n], h[n], pd[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_float_vec(v[i+2], w + i, h + i, pd + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h || pd[i] != 1)
			fail("bad %dth image size\n", i);
	}

	int out_pd = 8;
	float (*y) = xmalloc(*w * *h * out_pd * sizeof*y);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < *w * *h; i++)
	{
		float tmp[n];
		int ngood = 0;
		for (int j = 0; j < n; j++)
			if (isgood(x[j]+i, 1)) {
				tmp[ngood] = x[j][i];
				ngood += 1;
			}
		float_varkmeans(y + i*out_pd, out_pd, tmp, ngood, PRECISION());
		if (VECOH_VERBOSE() && 2*i == *w**h)
		{ // if requested, print the computation at the central pixel
			float *Y = y + i*out_pd;
			fprintf(stderr, "x = ");
			for (int l = 0; l < n; l++)
				fprintf(stderr, "%g%c",tmp[l],l==n-1?'\n':' ');
			fprintf(stderr, "\tn = %g", Y[0]);
			fprintf(stderr, "\ty = ");
			for (int l = 0; l < Y[0]; l++)
				fprintf(stderr,"%g%c",Y[1+l],*Y-1==l?'\n':' ');

		}
	}
	iio_write_image_float_vec(filename_out, y, *w, *h, out_pd);
	free(y);
	for (int i = 0; i < n; i++)
		free(x[i]);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_vecoh(c, v); }
#endif
