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

// y[k] = sum_i x[i][k]
static void float_sum(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k];
	}
}

// y[k] = (1/n) * sum_i x[i][k]
static void float_avg(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

// y[k] = prod_i x[i][k]
static void float_mul(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 1;
		for (int i = 0; i < n; i++)
			y[k] *= x[i][k];
	}
}

static float fnorm(float *x, int n)
{
	switch(n) {
	case 1: return fabs(x[0]);
	case 2: return hypot(x[0], x[1]);
	default: return hypot(x[0], fnorm(x+1, n-1));
	}
}

// y[] = smallest x[i][] in euclidean norm
static void float_min(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float mino = fnorm(x[0], d);
	for (int i = 1; i < n; i++)
	{
		float ni = fnorm(x[i], d);
		if (ni < mino) {
			midx = i;
			mino = ni;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];
}

// y[] = largest x[i][] in euclidean norm
static void float_max(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float mino = fnorm(x[0], d);
	for (int i = 1; i < n; i++)
	{
		float ni = fnorm(x[i], d);
		if (ni > mino) {
			midx = i;
			mino = ni;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];
}

static float medscore(float *xx, int idx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	float r = 0;
	for (int i = 0; i < n; i++)
		if (i != idx)
		{
			float v[d];
			for (int j = 0; j < d; j++)
				v[j] = x[idx][j] - x[i][j];
			r += fnorm(v, d);
		}
	return r;
}

// y[] = x[i][] which is closest to the euclidean median
static void float_med(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float misc = medscore(xx, 0, d, n);
	for (int i = 1; i < n; i++)
	{
		float si = medscore(xx, i, d, n);
		if (si < misc) {
			midx = i;
			misc = si;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];

}

static float float_mod_1d(float *x, int n)
{
	float h[0x100];
	for (int i = 0; i < 0x100; i++)
		h[i] = 0;
	for (int i = 0; i < n; i++)
	{
		int xi = x[i];
		if (xi < 0) fail("negative xi=%g", x[i]);//xi = 0;
		if (xi > 0xff) fail("large xi=%g", x[i]);//xi = 0xff;
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


// y[k] = mode of all x[i][k]
static void float_modc(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int i = 0; i < d; i++) {
		float t[n];
		for (int j = 0; j < n; j++)
			t[j] = x[j][i];
		y[i] = float_mod_1d(t, n);
	}
}


// euclidean distance between the vectors x and y
static float fdist(float *x, float *y, int n)
{
	return n ? hypot(*x - *y, fdist(x + 1, y + 1, n - 1)) : 0;
}

// euclidean distance between the vectors x and y, regularized around 0
static float fdiste(float *x, float *y, int n, float e)
{
	return n ? hypot(*x - *y, fdiste(x + 1, y + 1, n - 1, e)) : e;
}

#include "smapa.h"


int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

//static float float_med(float *x, int n)
//{
//	if (!n) return NAN;//fail("empty list of pixel values!");
//	if (n == 1) return x[0];
//	if (n == 2) return x[0];
//	qsort(x, n, sizeof*x, compare_floats);
//	return x[n/2];
//}

static float average_of_cluster_variances(
		int *group, float *mean, int k,
		float *x, int n)
{
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

	// run a few iterations of k-means
	for (int cx = 0; cx < 5; cx++)
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

SMART_PARAMETER_SILENT(VECOH_K,2)

// y[0] = number of modes of the samples x[i] (must satisfy y[0] < dy-1)
// y[1] = position of first mode
// y[2] = position of second mode
// ...
// y[y[0]] = position of the last mode
// y[y[0]+1] = nan
// ...
// y[dy-1] = nan
static void float_kmeans(float *y, int dy, float *x, int n)
{
	int k = VECOH_K();
	if (k >= dy) k = dy - 1;
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

static int main_kmeans(int c, char *v[])
{
	if(c<3)return fprintf(stderr,"usage:\n\tVECOH_K=2 kmeans x1 ... xn\n");
	int n = c - 1;
	float x[n], y[n];
	for (int i = 0; i < n; i++)
		x[i] = atof(v[i+1]);
	float_kmeans(y, 20, x, n);
	int k = y[0];
	//for (int i = 0; i < k; i++)
	//	fprintf(stderr, "center[%d] = %g\n", i+1, y[1+i]);
	int group[n];
	assign_each_point_to_nearest_mean(group, y+1, k, x, n);
	//for(int i = 0; i < n; i++)fprintf(stderr,"x[%d] = %g\n",i,x[i]);
	//for(int i = 0; i < n; i++)fprintf(stderr,"g[%d] = %d\n",i,group[i]);
	float acv = average_of_cluster_variances(group, y+1, k, x, n);
	fprintf(stderr, "racv = %g\t", sqrt(acv));
	fprintf(stderr, "got %d centers (", k);
	for (int i = 0; i < k; i++)
		fprintf(stderr, "%g%c", y[1+i], i+1==k?')':' ');
	puts("");
	return 0;
}

int main(int c, char **v) { return main_kmeans(c, v); }

static bool isgood(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isfinite(x[i]))
			return false;
	return true;
}

#include "pickopt.c"

//int main_vecov(int c, char *v[])
//{
//	char *filename_out = pick_option(&c, &v, "o", "-");
//	if (c < 4) {
//		fprintf(stderr,
//		"usage:\n\t%s {kmeans|modes} [v1 ...] [-o out]\n", *v);
//		//          0  1              2  3
//		return EXIT_FAILURE;
//	}
//	int n = c - 2;
//	char *operation_name = v[1];
//	void (*f)(float*,float*,int,int) = NULL;
//	if (0 == strcmp(operation_name, "kmeans")) f = float_kmeans;
//	if (0 == strcmp(operation_name, "modes"))  f = float_modes;
//	if (!f) fail("unrecognized operation \"%s\"", operation_name);
//	float *x[n];
//	int w[n], h[n], pd[n];
//	for (int i = 0; i < n; i++)
//		x[i] = iio_read_image_float_vec(v[i+2], w + i, h + i, pd + i);
//	for (int i = 0; i < n; i++) {
//		if (w[i] != *w || h[i] != *h || pd[i] != *pd)
//			fail("%dth image sizes mismatch\n", i);
//	}
//	int out_pd = *pd;
//	float (*y) = xmalloc(*w * *h * out_pd * sizeof*y);
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//	for (int i = 0; i < *w * *h; i++) {
//		float tmp[n][*pd];
//		int ngood = 0;
//		for (int j = 0; j < n; j++)
//			if (isgood(x[j]+i**pd, *pd)) {
//				for (int k = 0; k < *pd; k++)
//					tmp[ngood][k] = x[j][i**pd+k];
//				ngood += 1;
//			}
//		f(y + i**pd, tmp[0], *pd, ngood);
//	}
//	iio_write_image_float_vec(filename_out, y, *w, *h, *pd);
//	free(y);
//	for (int i = 0; i < n; i++)
//		free(x[i]);
//	return EXIT_SUCCESS;
//}

//#ifndef HIDE_ALL_MAINS
//int main(int c, char **v) { return main_vecov(c, v); }
//#endif
