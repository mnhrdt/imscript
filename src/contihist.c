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

void fill_continuous_histogram_fake(long double (*bins)[2], int nbins,
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

int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

// thi
static void sort_quad(struct quad *q)
{
	qsort(q->abcd, 4, sizeof q->abcd[0], compare_floats);
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

static int degenerateP(float *x)
{
	return !(x[0] < x[1] && x[1] < x[2] && x[2] < x[3]);
}

void fill_continuous_histogram(long double (*bins)[2], int nbins,
		float min, float max, float *x, int w, int h)
{
	// initialize bins
	for (int i = 0; i < nbins; i++)
	{
		float idi = min + i * (max - min) / nbins;
		bins[i][0] = idi;
		bins[i][1] = 0;
	}

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
		if (degenerateP(abcd)) continue;
		nquads += 1;

		// compute bins associated to this quad
		int iabcd[4];
		for (int k = 0; k < 4; k++)
			iabcd[k] = ifloatbin(nbins, min, max, abcd[k]);

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
			e[iabcd[k]] += dslope_abcd[k];
	}

	fprintf(stderr, "generic quads = %d / %d (%g%%)\n",
			nquads, (w-1)*(h-1), nquads*100.0/((w-1)*(h-1)));

	// integrate slope increments (twice)
	for (int i = 1; i < nbins; i++) e[i] += e[i-1];
	e[0] = 0;
	for (int i = 1; i < nbins; i++) e[i] += e[i-1];
	e[0] = 0;


	// copy to output
	for (int i = 0; i < nbins; i++)
		bins[i][1] = e[i];


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

void accumulate_histogram(long double (*h)[2], int n)
{
	for (int i = 1; i < n; i++)
		h[i][1] += h[i-1][1];
}

void dump_histogram(long double (*h)[2], int n)
{
	long double (*a)[2] = xmalloc(n*sizeof*a);
	memcpy(a, h, n*sizeof*a);
	accumulate_histogram(a, n);
	printf("set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
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
		printf("unset key\n");
	printf("plot \"-\" w impulses title \"histogram\", \"-\" w lines title \"accumulated histogram\",0");
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
	long double f = m/a[n-1][1];
	//fprintf(stderr, "m = %Lg\n", m);
	//fprintf(stderr, "f = %Lg\n", f);
	//fprintf(stderr, "a = %Lg\n", a[n-1][1]);
	//fprintf(stderr, "n = %d\n", n);
	for (int i = 0; i < n; i++)
		printf("\t%Lg\t%Lg\n", a[i][0], f*a[i][1]);
	printf("end\n");
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


#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s nbins min max img.png\n", *v);
		//                          0 1     2   3   4
	}
	int nbins = atoi(v[1]);
	float hmin = atof(v[2]);
	float hmax = atof(v[3]);
	char *filename_in = v[4];

	// read input image
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	// allocate space for the histogram data
	long double bins[nbins][2];

	// compute continuous histogram
	fill_continuous_histogram(bins, nbins, hmin, hmax, x, w, h);

	// dump histogram to stdout
	dump_histogram(bins, nbins);

	// cleanup and exit
	free(x);
	return 0;
}
