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
#include "modes_detector.c"

#include "smapa.h"
SMART_PARAMETER(PRECISION,10)
SMART_PARAMETER_SILENT(NUMBER_OF_KMEANS_ITERATIONS,5)


static int compare_floats(const void *a, const void *b)
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

// compute the avg of the diq of each cluster
static float average_of_cluster_diq(
		int *group, float *mean, int k,
		float *x, int n)
{
	mean;
	if (!k) return 0;
	float q25[k], q75[k];

	for (int j = 0; j < k; j++)
	{
		int tn = 0;
		float tx[n];
		for (int i = 0; i < n; i++)
			if (j == group[i])
				tx[tn++] = x[i];
		assert(tn);
		int a[2] = {(tn+1)/4-1, tn/4 };
		int b[2] = {(3*tn+1)/4-1, (3*tn)/4 };
		//fprintf(stderr, "\tq25(%d) = %d %d\tq75(%d) = %d %d\n",
		//		tn, a[0], a[1], tn, b[0], b[1]);
		for (int i = 0; i < 2; i++) {
			if (a[0] < 0) a[0] = 1;
			if (a[1] < 0) a[1] = 1;
			if (b[0] < 0) b[0] = 1;
			if (b[1] < 0) b[1] = 1;
			if (a[0] >= tn) a[0] = tn-1;
			if (a[1] >= tn) a[1] = tn-1;
			if (b[0] >= tn) b[0] = tn-1;
			if (b[1] >= tn) b[1] = tn-1;
			assert(a[i] >= 0); assert(b[i] >= 0);
			assert(a[i] < tn); assert(b[i] < tn);
		}
		// TODO: use the corrected coefficients of 1/3 & 2/3 here:
		q25[j] = (tx[a[0]] + tx[a[1]])/2;
		q75[j] = (tx[b[0]] + tx[b[1]])/2;
	}

	float r = 0;
	for (int i = 0; i < k; i++)
		r += q75[i] - q25[i];
	return r / k;
}

// given a set of centers, assign the index of the nearest center to each point
static void assign_each_point_to_nearest_center(int *assign, float *center,
		int k, float *x, int n)
{
	// TODO: this function should run in O(n) instead of O(n*k)
	// using the fact than the x[i] are ordered
	// and the assignement is thus monotonic
	for (int i = 0; i < n; i++)
	{
		int bj = 0;
		for (int j = 0; j < k; j++)
			if ( fabs(x[i] - center[j]) < fabs(x[i] - center[bj]) )
				bj = j;
		assign[i] = bj;
	}
}

// INPUT  x[n] : ordered array of floats
// INPUT  g[n] : array of (0..k), saying to which cluster belongs each number
// OUTPUT m[k] : ordered array of cluster centers
static void update_means_of_each_group(float *m, int *g, int k, float *x, int n)
{
	int c[k];
	for (int j = 0; j < k; j++)
		m[j] = c[j] = 0;
	for (int i = 0; i < n; i++)
	{
		m[g[i]] += x[i];
		c[g[i]] += 1;
	}
	for (int j = 0; j < k; j++)
		if (c[j])
			m[j] /= c[j];
		else
			m[j] = NAN;
}

static
void update_medians_of_each_group(float *m, int *g, int k, float *x, int n)
{
	for (int i = 0; i < n-1; i++)
		assert(x[i] <= x[i+1]);
	//fprintf(stderr, "\tupdating medians %d%d :\n", k, n);
	//fprintf(stderr, "\tx : ");
	//for (int i = 0; i < n; i++)
	//	fprintf(stderr, "%g%c", x[i], i==n-1?'\n':' ');
	//fprintf(stderr, "\tg : ");
	//for (int i = 0; i < n; i++)
	//	fprintf(stderr, "%d%c", g[i], i==n-1?'\n':' ');

	for (int j = 0; j < k; j++)
	{
		int tn = 0;
		float tx[n];
		for (int i = 0; i < n; i++)
			if (j == g[i])
				tx[tn++] = x[i];
		if (tn) {
			int a = (tn+1)/2-1;
			int b = tn/2;
			assert(a >= 0); assert(b >= 0);
			assert(a < tn); assert(b < tn);
			m[j] = (tx[a] + tx[b])/2;
			//fprintf(stderr, "\t\tnt a b m[%d] = %d %d %d %g\n",
			//		j, tn, a, b, m[j]);
		} else
			m[j] = NAN;
	}
}

// This function sorts and processes a vector x[n] for the trivial cases of 1d
// clustering algorithms.  It returns 1 if no further processing is needed.
//
// x[n] : array of floats that will be sorted
// m[k] : computed clusters or nan, in some trivial cases
// return value : whether the compoutation is already finished for this vector.
static int sort_and_sanitize_and_check(float *m, float *x, int n, int k)
{
	for (int i = 0; i < k; i++)
		m[i] = NAN;
	if (!n)
		return 1;
	qsort(x, n, sizeof*x, compare_floats);
	if (n <= k) {
		for (int i = 0; i < n; i++)
			m[i] = x[i];
		return 1;
	}
	return 0;
}

// INPUT  x[n] : ordered array of numbers
// OUTPUT m[k] : centers uniformly-distributed clusters, to be filled-in
static void initialize_clusters_uniformly(float *m, float *x, int n, int k)
{
	assert(n > k);
	for (int i = 0; i < k; i++)
	{
		int idx = lrint(  (i + 1.0) * (n - 1.0) / (k + 1.0)  );
		assert(idx >= 0);
		assert(idx < n);
		m[i] = x[idx];
	}

	// The following hack is needed for the particular case when
	// a lot of the input numbers are exactly the same number.
	//
	// It solves the case when all the clusters collapse into this point.
	float wiggle = 0.00001;
	for (int i = 1; i < k; i++)
		if (m[i-1] == m[i])
		{
			m[i-1] -= wiggle;
			m[i] += wiggle;
		}
}

// the standard k-means algorithm for 1d points
static void kmeans_1d(float *means, float *x, int n, int k)
{
	if (sort_and_sanitize_and_check(means, x, n, k))
		return;

	initialize_clusters_uniformly(means, x, n, k);

	int numit = NUMBER_OF_KMEANS_ITERATIONS();
	for (int cx = 0; cx < numit; cx++)
	{
		int group[n]; // to which cluster belongs each point
		assign_each_point_to_nearest_center(group, means, k, x, n);
		update_means_of_each_group(means, group, k, x, n);
	}
}

// the standard k-means algorithm for 1d points
static void kmedians_1d(float *means, float *x, int n, int k)
{
	if (sort_and_sanitize_and_check(means, x, n, k))
		return;

	initialize_clusters_uniformly(means, x, n, k);

	int numit = NUMBER_OF_KMEANS_ITERATIONS();
	for (int cx = 0; cx < numit; cx++)
	{
		int group[n]; // to which cluster belongs each point
		assign_each_point_to_nearest_center(group, means, k, x, n);
		update_medians_of_each_group(means, group, k, x, n);
		//update_means_of_each_group(means, group, k, x, n);
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

static void float_kmedians(float *y, int dim_y, float *x, int n, int k)
{
	for (int i = 0; i < dim_y; i++) y[i] = NAN;
	if (k >= dim_y) k = dim_y - 1;
	float medians[k];
	kmedians_1d(medians, x, n, k);
	int count = 0;
	for (int i = 0; i < k; i++)
		count += isfinite(medians[i]);
	y[0] = count;
	int ypos = 1;
	for (int i = 0; i < k; i++)
		if (isfinite(medians[i]))
			y[ypos++] = medians[i];
}

// run k-means with increasing k until the average variance of the clusters
// falls below a pre-sed precision
static void float_xmeans(float *y, int dim_y, float *x, int n, float prec)
{
	for (int pre_k = 1; pre_k < dim_y; pre_k++)
	{
		float_kmeans(y, dim_y, x, n, pre_k);
		int k = y[0];
		int group[n];
		assign_each_point_to_nearest_center(group, y+1, k, x, n);
		float acv = average_of_cluster_variances(group, y+1, k, x, n);
		if (acv < prec*prec)
			break;
	}
}

// run k-means with increasing k until the average variance of the clusters
// falls below a pre-sed precision
static void float_xmedians(float *y, int dim_y, float *x, int n, float prec)
{
	for (int pre_k = 1; pre_k < dim_y; pre_k++)
	{
		//fprintf(stderr,"xmedians of n = %d with pre_k = %d\n",n,pre_k);
		float_kmedians(y, dim_y, x, n, pre_k);
		int k = y[0];
		int group[n];
		assign_each_point_to_nearest_center(group, y+1, k, x, n);
		float acv = average_of_cluster_diq(group, y+1, k, x, n);
		//float acv = average_of_cluster_variances(group, y+1, k, x, n);
		//fprintf(stderr, "...got k=%d acv=%g\n", k, acv);
		//for (int i = 0; i < n; i++)
		//	fprintf(stderr, "...x[%d] = %g\n", i, x[i]);
		//for (int i = 0; i < n; i++)
		//	fprintf(stderr, "...group[%d] = %d\n", i, group[i]);
		//for (int i = 0; i < 1+k; i++)
		//	fprintf(stderr, "...y[%d] = %g\n", i, y[i]);
		if (acv < prec)
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
	assign_each_point_to_nearest_center(group, y+1, k, x, n);
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
	if(c<4)return fprintf(stderr,"usage:\n\t"
			"vecoh {kmeans|kmedians} x1 ... xn\n");

	char * operation_name = v[1];
	void (*f)(float*,int,float*,int,float) = NULL;
	if (0==strcmp(operation_name,"kmeans"))    f=float_xmeans;
	if (0==strcmp(operation_name,"kmedians"))  f=float_xmedians;
	if (0==strcmp(operation_name,"contrario")) f=acontrario_modes_detector;
	if (!f) return fprintf(stderr, "xunrecognized operation \"%s\"\n",
				operation_name);

	int n = c - 2;
	int max_k = 20;
	float x[n], y[max_k];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[i+2]);
	f(y, max_k, x, n, PRECISION());
	int k = y[0];
	int group[n];
	assign_each_point_to_nearest_center(group, y+1, k, x, n);
	if (f == float_xmeans) {
		float acv = average_of_cluster_variances(group, y+1, k, x, n);
		fprintf(stderr, "racv = %g\t", sqrt(acv));
	}
	if (f == float_xmedians) {
		float acv = average_of_cluster_diq(group, y+1, k, x, n);
		fprintf(stderr, "racv = %g\t", acv);
	}
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
	if (c < 3) {
		fprintf(stderr,
		"usage:\n\t%s {kmeans|kmedians|contrario} [v1 ...] [-o out]\n", *v);
		//          0  1                            2  3
		return 1;
	}
	int n = c - 2;
	char *operation_name = v[1];
	float *x[n];
	int w[n], h[n], pd[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_float_vec(v[i+2], w + i, h + i, pd + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h || pd[i] != 1)
			fail("bad %dth image size\n", i);
	}

	void (*f)(float*,int,float*,int,float) = NULL;
	if (0==strcmp(operation_name, "kmeans"))    f=float_xmeans;
	if (0==strcmp(operation_name, "kmedians"))  f=float_xmedians;
	if (0==strcmp(operation_name, "contrario")) f=acontrario_modes_detector;
	if (!f) return fprintf(stderr, "unrecognized operation \"%s\"\n",
				operation_name);

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
		f(y + i*out_pd, out_pd, tmp, ngood, PRECISION());
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
