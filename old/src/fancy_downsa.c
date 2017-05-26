#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "fancy_image.h"

static double average_of_non_nans(double *x, int n)
{
	int cx = 0;
	double ax = 0;
	for (int i = 0; i < n; i++)
	{
		if (isnan(x[i])) continue;
		cx += 1;
		ax += x[i];
	}
	return cx ? ax / cx : NAN;
}

static double combine_4doubles(double v[4], int m)
{
	if (m == 'f') return v[0];
	if (m == 'v') return (v[0]+v[1]+v[2]+v[3])/4;
	if (m == 'V') return average_of_non_nans(v, 4);
	if (m == 'i') return fmin(fmin(v[0],v[1]), fmin(v[2],v[3]));
	if (m == 'a') return fmax(fmax(v[0],v[1]), fmax(v[2],v[3]));
	return m;//NAN;
}

static float count_non_nans(float *v, int n)
{
	float r = 0;
	for (int i = 0; i < n; i++)
		r += !isnan(v[i]);
	return r;
}

static float min_non_nans(float *v, int n)
{
	float r = INFINITY;
	for (int i = 0; i < n; i++)
		if (!isnan(v[i]))
			r = fmin(r, v[i]);
	return r;
}

static float max_non_nans(float *v, int n)
{
	float r = -INFINITY;
	for (int i = 0; i < n; i++)
		if (!isnan(v[i]))
			r = fmax(r, v[i]);
	return r;
}

static float avg_non_nans(float *v, int n)
{
	float r = 0;
	int nr = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(v[i]))
		{
			r += v[i];
			nr += 1;
		}
	return nr ? r / nr : NAN;
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static float med_non_nans(float *v, int n)
{
	float t[n];
	for (int i = 0; i < n; i++)
		t[i] = v[i];
	qsort(t, n, sizeof*t, compare_floats);
	return 0.5 * (t[n/2-1] + t[n/2]);
}

static float combine_floats(float *v, int n, int op)
{
	switch(op) {
	case 'f': return *v;
	case 'k': return count_non_nans(v, n);
	case 'i': return min_non_nans(v, n);
	case 'a': return max_non_nans(v, n);
	case 'v': return avg_non_nans(v, n);
	case 'e': return med_non_nans(v, n);
	default: return NAN;
	}
}

static void fancy_downsa(char *fname_out, char *fname_in, int nw, int nh, int m)
{
	// open input image
	struct fancy_image *a = fancy_image_open(fname_in, "r");

	// read information from input image
	int tw = 0, th = 0, fmt = 0, bps = 0;
	int tiffo = fancy_image_leak_tiff_info(&tw, &th, &fmt, &bps, a);

	// create output image of the appropriate size and options
	int bw = a->w / nw;
	int bh = a->h / nh;
	struct fancy_image *b = fancy_image_create(fname_out,
			"w=%d,h=%d,pd=%d,bps=%d,fmt=%d,tw=%d,th=%d",
			bw, bh, a->pd, bps, fmt, tw, th);

	// fill-in the zoomed-out image
	for (int j = 0; j < b->h; j++)
	for (int i = 0; i < b->w; i++)
	for (int l = 0; l < b->pd; l++)
	{
		int nv = 0;
		float vv[nw*nh];
		for (int dj = 0; dj < nh; dj++)
		for (int di = 0; di < nw; di++)
		{
			int ii = i * nw + di;
			int jj = j * nh + dj;
			vv[nv++] = fancy_image_getsample(a, ii, jj, l);
		}
		float r = combine_floats(vv, nv, m);
		fancy_image_setsample(b, i, j, l, r);
	}

	// close both images
	fancy_image_close(b);
	fancy_image_close(a);
}

#include <stdio.h>
int main_fancy_downsa(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t"
			"%s {i|e|a|v|r|f} pw ph in.tiff out.tiff\n", *v);
		//        0  1            2  3  4       5
		return 1;
	}
	int op = v[1][0];
	int pw = atoi(v[2]);
	int ph = atoi(v[3]);
	char *filename_in = v[4];
	char *filename_out = v[5];

	fancy_downsa(filename_out, filename_in, pw, ph, op);

	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_fancy_downsa(c, v); }
#endif
