#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static float floatbin(int nbins, float min, float max, float x)
{
	if (isnan(x)) return x;
	if (x < min) x = min;
	if (x > max) x = max;
	float f = (nbins-1) * (x - min) / (max - min);
	return f;
}

static int ifloatbin(int nbins, float min, float max, float x)
{
	assert(isfinite(x));
	if (x < min) x = min;
	if (x > max) x = max;
	float f = (nbins-1) * (x - min) / (max - min);
	int r = lrint(f);
	assert(r >= 0);
	assert(r < nbins);
	return r;
}


#include "xfopen.c"
static void odump(char *fname, long double (*bin)[2], int n)
{
	FILE *f = xfopen(fname, "w");
	for (int i = 0; i < n; i++)
		fprintf(f, "%d %Lf %Lf\n", i, bin[i][0], bin[i][1]);
	xfclose(f);
}

static void edump(char *fname, long double *e, int n)
{
	FILE *f = xfopen(fname, "w");
	for (int i = 0; i < n; i++)
		fprintf(f, "%d %Lf\n", i, e[i]);
	xfclose(f);
}

static void fill_continuous_histogram_fake(long double (*bins)[2], int nbins,
		float min, float max, float *x, int w, int h)
{
	// initialize bins
	for (int i = 0; i < nbins; i++)
	{
		float idi = min + i * (max - min) / nbins;
		bins[i][0] = idi;
		bins[i][1] = 0;
	}

	// accumulate
	for (int i = 0; i < w*h; i++)
	{
		float fi = floatbin(nbins, min, max, x[i]);
		if (isfinite(fi))
		{
			assert(fi >= 0);
			assert(fi < nbins);
			int ifi = floor(fi);
			assert(ifi >= 0);
			assert(ifi < nbins);
			bins[ifi][1] += 1;
		}
	}
}

#include "xmalloc.c"

// a quad is a square cell bounded by 4 pixels
// it knows its four values and its position on the original image
struct quad {
	float abcd[4];
	int ij;
};

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}


static void sort_quad(struct quad *q)
{
	qsort(q->abcd, 4, sizeof q->abcd[0], compare_floats);
}

static void sort_four_values(float *x)
{
	qsort(x, 4, sizeof*x, compare_floats);
}


//// an event indicates that the histogram changes slope at e->q->abcd[e->type]
//struct event {
//	struct quad *q;
//	float A;
//	int type;
//	int bin;
//};
//
//int compare_events(const void *a, const void *b)
//{
//	const struct event *ea = (const struct event *) a;
//	const struct event *eb = (const struct event *) b;
//	float fa = ea->q->abcd[ea->type];
//	float fb = eb->q->abcd[eb->type];
//	return (fa > fb) - (fa < fb);
//}

#include <stdbool.h>
static int degenerateP(float min, float max, float *x)
{
	if (x[0] < min) return true;
	if (x[0] > max) return true;
	return !(x[0] < x[1] && x[1] < x[2] && x[2] < x[3]);
}

void fill_continuous_histogram(long double (*bins)[2], int nbins,
		float min, float max, float *x, int w, int h)
{
	// initialize bins
	for (int i = 0; i < nbins; i++)
	{
		long double idi = min + i * (max - min) / nbins;
		bins[i][0] = idi;
		bins[i][1] = 0;
	}
	odump("/tmp/dch_full_0.txt", bins, nbins);

	// fill quads
	struct quad *q = xmalloc(w * h * sizeof*q);
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		int ij = j*w + i;
		q[ij].ij = ij;
		q[ij].abcd[0] = x[ij];
		q[ij].abcd[1] = x[ij+1];
		q[ij].abcd[2] = x[ij+w];
		q[ij].abcd[3] = x[ij+1+w];
		sort_quad(q + ij);
	}

	// fill events (each quad determines four events)
	long double *e  = xmalloc(nbins * sizeof*e);
	for (int i = 0; i < nbins; i++)
		e[i] = 0;
	int nquads = 0;
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		// discard degenerate quads
		float *abcd = q[j*w + i].abcd;
		if (degenerateP(min, max, abcd)) continue;
		nquads += 1;

		// compute bins associated to this quad
		int iabcd[4];
		for (int k = 0; k < 4; k++)
			iabcd[k] = ifloatbin(nbins, min, max, abcd[k]);

		//fprintf(stderr, "full cell (%g %g %g %g)[%d %d %d %d]\n",
		//		abcd[0], abcd[1], abcd[2], abcd[3],
		//		iabcd[0], iabcd[1], iabcd[2], iabcd[3]);

		// compute parameters associated to this quad
		float gamma = 2 / (abcd[3] + abcd[2] - abcd[1] - abcd[0]);
		float alpha = gamma / (abcd[1] - abcd[0]);
		float beta  = gamma / (abcd[2] - abcd[3]);
		assert(gamma > 0);
		assert(alpha > 0);
		assert(beta < 0);

		// accumulate slope increments
		float dslope_abcd[4] = { alpha, -alpha, beta, -beta};
		for (int k = 0; k < 4; k++)
		{
			e[iabcd[k]] += dslope_abcd[k];
		}
	}
	edump("/tmp/dch_full_1.txt", e, nbins);

	fprintf(stderr, "generic quads = %d / %d (%g%%)\n",
			nquads, (w-1)*(h-1), nquads*100.0/((w-1)*(h-1)));

	// integrate slope increments (twice)
	for (int i = 1; i < nbins; i++) e[i] += e[i-1];
	e[0] = 0;
	edump("/tmp/dch_full_2.txt", e, nbins);
	for (int i = 1; i < nbins; i++) e[i] += e[i-1];
	e[0] = 0;
	edump("/tmp/dch_full_3.txt", e, nbins);


	// copy to output
	for (int i = 0; i < nbins; i++)
		bins[i][1] = e[i];
	odump("/tmp/dch_full_4.txt", bins, nbins);


//	struct event *e = xmalloc(w * h * 4 * sizeof*e);
//	int ecx = 0;
//	for (int j = 0; j < h - 1; j++)
//	for (int i = 0; i < w - 1; i++)
//	{
//		int ij = j*w + i;
//		for (int k = 0; k < 4; k++)
//		{
//			float *abcd = e[ecx].q->abcd;
//			e[ecx].q = q + ij;
//			e[ecx].type = k;
//			// TODO: copy appropriate formula here
//			e[ecx].A = 0;
//			ecx += 1;
//		}
//	}

	// construct histogram
	// (the events describe the derivative of the histogram, thus to
	// compute the histogram, we have to integrate the events)

	//// sort events
	//qsort(e, ecx, sizeof*e, compare_events);

	//// construct histogram
	//float A = 0, B = 0;
	//int eidx = 0;
	//for (int i = 0; i < nbins - 1; i++)
	//{
	//	// accumulate all the events that happen on this bin
	//	while (eidx < ecx)
	//	{
	//		struct event *ev = e + eidx;
	//		if (ev->q->abcd[ev->type] >= bins[i+1][0])
	//			break;
	//		A += ev->A;
	//		B += ev->B;
	//		eidx += 1;
	//	}
	//	bins[i][1] = A * bins[i][0] + B;
	//}

}

// copy the values of a square cell into an array, for easy access
static void copy_cell_values(float q[4], float *x, int w, int h, int i, int j)
{
	assert(0 <= i); assert(i < w - 1);
	assert(0 <= j); assert(j < h - 1);
	q[0] = x[ (j+0)*w + (i+0) ];
	q[1] = x[ (j+0)*w + (i+1) ];
	q[2] = x[ (j+1)*w + (i+0) ];
	q[3] = x[ (j+1)*w + (i+1) ];
}

// obtain the histogram bin that corresponds to the given value
static int bin(int n, float m, float M, float x)
{
	float f = (n - 1) * (x - m) / (M - m);
	int r = lrint(f);
	assert(m <= x); assert(x <= M); assert(0 <= r); assert(r < n);
	return r;
}

static int cell_is_degenerate(float m, float M, float q[4])
{
	for (int i = 0; i < 4; i++)
		if (q[i] < m || q[i] > M)
			return 1;
	return !(q[0] < q[1] && q[1] < q[2] && q[2] < q[3]);
}

static void integrate_values(long double (*o)[2], int n)
{
	// TODO : multiply each increment by the span of the interval
	for (int i = 1; i < n; i++)
		o[i][1] += o[i-1][1];// * (o[i][0] - o[i-1][0]);
	//o[0][1] = 0;
}

static void accumulate_jumps_for_one_cell(long double (*o)[2],
		int n, float m, float M, float q[4])
{
	// discard degenerate cells
	if (cell_is_degenerate(m, M, q))
		return;

	// give nice names to numbers
	float A = q[0];
	float B = q[1];
	float C = q[2];
	float D = q[3];
	int i_A = bin(n, m, M, A);
	int i_B = bin(n, m, M, B);
	int i_C = bin(n, m, M, C);
	int i_D = bin(n, m, M, D);
	A = i_A;
	B = i_B;
	C = i_C;
	D = i_D;

	if (i_A == i_B || i_B == i_C || i_C == i_D) return;
	//fprintf(stderr, "simp cell (%g %g %g %g)[%d %d %d %d]\n",
	//		A, B, C, D, i_A, i_B, i_C, i_D);

	assert(i_A < i_B);
	assert(i_B < i_C);
	assert(i_C < i_D);

	double fac = 1;//33333;

	// accumulate jumps
	o[ i_A ][1] += fac * 2 / (C + D - B - A) / (B - A);
	o[ i_B ][1] -= fac * 2 / (C + D - B - A) / (B - A);
	o[ i_C ][1] -= fac * 2 / (C + D - B - A) / (D - C);
	o[ i_D ][1] += fac * 2 / (C + D - B - A) / (D - C);

	//fprintf(stderr, "\tacc %d %lf\n", i_A, + fac*2/(C + D - B - A)/(B - A));
	//fprintf(stderr, "\tacc %d %lf\n", i_B, - fac*2/(C + D - B - A)/(B - A));
	//fprintf(stderr, "\tacc %d %lf\n", i_C, - fac*2/(C + D - B - A)/(D - C));
	//fprintf(stderr, "\tacc %d %lf\n", i_D, + fac*2/(C + D - B - A)/(D - C));
}


void fill_continuous_histogram_simple(
	long double (*o)[2], // output histogram array of (value,density) pairs
	int n,               // requested number of bins for the histogram
	float m,             // requested minimum of the histogram
	float M,             // requested maximum of the histogram
	float *x,            // input image data
	int w,               // input image width
	int h                // input image height
	)
{
	// initialize bins (n equidistant points between m and M)
	for (int i = 0; i < n; i++)
	{
		o[i][0] = m + i * (M - m) / n;
		o[i][1] = 0;
	}

	// compute 2nd derivative of histogram
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		float q[4]; // here we store the 4 values of the cell at (i,j)
		copy_cell_values(q, x, w, h, i, j);
		sort_four_values(q);
		accumulate_jumps_for_one_cell(o, n, m, M, q);
	}

	// integrate twice
	integrate_values(o, n);
	integrate_values(o, n);
}

void accumulate_histogram(long double (*h)[2], int n)
{
	for (int i = 1; i < n; i++)
		h[i][1] += h[i-1][1];
}

void dump_histogram(FILE *f, long double (*h)[2], int n)
{
	long double (*a)[2] = xmalloc(n*sizeof*a);
	memcpy(a, h, n*sizeof*a);
	accumulate_histogram(a, n);
	fprintf(f, "set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
	//printf("set yrange [0:]\n");
	//printf("set format y \"\"\n");
	//if (SHOWSTATS() > 0) {
	//	printf("set samples 1000\nset key left\n");
	//	printf("N(m,s,x)=exp(-(x-m)**2/(2*s*s))/(s*sqrt(2*pi))\n");
	//	printf("L(m,s,x)=exp(-sqrt(2)*abs(x-m)/s)/(s*sqrt(2))\n");
	//	printf("U0(a,b,x)=x<a?0:(x>b?0:1/(b-a))\n");
	//	printf("U(m,s,x)=U0(m-s*sqrt(3),m+s*sqrt(3),x)\n");
	//	//printf("C(m,s,x)=1/((1+((x-m)/s)**2)*(pi*s))\n");
	//}
	//else
		fprintf(f, "unset key\n");
	fprintf(f, "plot \"-\" w impulses title \"histogram\", \"-\" w lines title \"accumulated histogram\",0");
	//if (SHOWSTATS() > 0) {
	//	void print_gnuplot_stats_string(long double (*)[2], int, float);
	//	print_gnuplot_stats_string(h, n, 1/SHOWSTATS());
	//}
	fprintf(f, "\n");
	long double m = 0;
	for (int i = 0; i < n; i++)
	{
		if (h[i][1] > m) m = h[i][1];
		fprintf(f, "\t%Lg\t%Lg\n", h[i][0], h[i][1]);
	}
	fprintf(f, "end\n");
	long double fac = m/a[n-1][1];
	//fprintf(stderr, "m = %Lg\n", m);
	//fprintf(stderr, "f = %Lg\n", f);
	//fprintf(stderr, "a = %Lg\n", a[n-1][1]);
	//fprintf(stderr, "n = %d\n", n);
	for (int i = 0; i < n; i++)
		fprintf(f, "\t%Lg\t%Lg\n", a[i][0], fac*a[i][1]);
	fprintf(f, "end\n");
	free(a);
}

void dump_histogram_noacc(long double (*h)[2], int n)
{
	//long double (*a)[2] = xmalloc(n*sizeof*a);
	//memcpy(a, h, n*sizeof*a);
	//accumulate_histogram(a, n);
	printf("set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
	printf("set yrange [0:]\n");
	printf("set format y \"\"\n");
	//if (SHOWSTATS() > 0) {
	//	printf("set samples 1000\nset key left\n");
	//	printf("N(m,s,x)=exp(-(x-m)**2/(2*s*s))/(s*sqrt(2*pi))\n");
	//	printf("L(m,s,x)=exp(-sqrt(2)*abs(x-m)/s)/(s*sqrt(2))\n");
	//	printf("U0(a,b,x)=x<a?0:(x>b?0:1/(b-a))\n");
	//	printf("U(m,s,x)=U0(m-s*sqrt(3),m+s*sqrt(3),x)\n");
	//	//printf("C(m,s,x)=1/((1+((x-m)/s)**2)*(pi*s))\n");
	//}
	//else
		printf("unset key\n");
	//printf("plot \"-\" w impulses title \"histogram\", \"-\" w lines title \"accumulated histogram\"");
	printf("plot \"-\" w impulses");
	//if (SHOWSTATS() > 0) {
	//	void print_gnuplot_stats_string(long double (*)[2], int, float);
	//	print_gnuplot_stats_string(h, n, 1/SHOWSTATS());
	//}
	printf("\n");
	long double m = 0;
	for (int i = 0; i < n; i++)
	{
		if (h[i][1] > m) m = h[i][1];
		printf("\t%Lg\t%Lg\n", h[i][0], h[i][1]);
	}
	printf("end\n");
	//long double f = m/a[n-1][1];
	//fprintf(stderr, "m = %Lg\n", m);
	//fprintf(stderr, "f = %Lg\n", f);
	//fprintf(stderr, "a = %Lg\n", a[n-1][1]);
	//fprintf(stderr, "n = %d\n", n);
	//for (int i = 0; i < n; i++)
	//	printf("\t%Lg\t%Lg\n", a[i][0], f*a[i][1]);
	//printf("end\n");
	//free(a);
}

static void update_min_max_if_not_finite(float *m, float *M, float *x, int n)
{
	if (isfinite(*m) && isfinite(*M)) return;

	float min = INFINITY;
	float max = -INFINITY;
	for (int i = 0; i < n; i++)
	{
		if (x[i] < min) min = x[i];
		if (x[i] > max) max = x[i];
	}
	if (!isfinite(*m)) *m = min;
	if (!isfinite(*M)) *M = max;
}


#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main_contihist(int c, char *v[])
{
	// process input arguments
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s nbins min max img.png\n", *v);
		//                          0 1     2   3   4
		return 1;
	}
	int nbins = atoi(v[1]);
	float hmin = atof(v[2]);
	float hmax = atof(v[3]);
	char *filename_in = v[4];

	// read input image
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	// allocate space for the histogram data
	long double bins[3+nbins][2];

	update_min_max_if_not_finite(&hmin, &hmax, x, w*h);

	// compute continuous histogram
	//fill_continuous_histogram(bins, nbins, hmin, hmax, x, w, h);
	fill_continuous_histogram_simple(bins, nbins, hmin, hmax, x, w, h);

	// dump histogram to stdout
	dump_histogram(stdout, bins, nbins);

	// cleanup and exit
	free(x);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_contihist(c, v); }
#endif
