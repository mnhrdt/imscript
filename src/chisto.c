#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"

#include "xmalloc.c"
#include "fail.c"
#include "smapa.h"

SMART_PARAMETER_SILENT(HISMOTTH,0)

static void smooth_histogram_rw(long double (*h)[2], int n, int w)
{
	for (int i = 0; i < n-w; i++)
	{
		long double a = 0;
		for (int j = 0; j <= w; j++)
			a += h[i+j][1];
		h[i][1] = a;
	}
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

int fill_histogram(long double (*h)[2], float *in_x, int n)
{
	if (!n) return 0;
	float *x = xmalloc(n * sizeof*x);
	for (int i = 0; i < n; i++)
		x[i] = in_x[i];
	qsort(x, n, sizeof*x, compare_floats);
	h[0][0] = x[0];
	h[0][1] = 1;
	int r = 0;
	for (int i = 1; i < n; i++) {
		if (x[i] != x[i-1]) {
			r += 1;
			h[r][0] = x[i];
			h[r][1] = 0;
		}
		h[r][1] += 1;
	}
	r += 1;

	free(x);
	if (HISMOTTH() > 0)
		smooth_histogram_rw(h, r, HISMOTTH());
	return r;
}

void accumulate_histogram(long double (*h)[2], int n)
{
	for (int i = 1; i < n; i++)
		h[i][1] += h[i-1][1];
}

void dump_histograms(long double (*h[4])[2], int n[4])
{
	printf("set terminal pngcairo size 800, 600\n");
	printf("set style fill transparent solid 0.33 border\n");
	printf("#set style fill transparent pattern 4 border\n");
	printf("set style function filledcurves y1=0\n");
	printf("set style data filledcurves y1=0\n");
	printf("set xrange [0:255]\n");
	printf("set yrange [0:]\n");
	printf("set format y \"\"\n");
	printf("unset key\n");
	printf("plot "
			"\"-\" lw 1 lc rgb \"red\", "
			"\"-\" lw 1 lc rgb \"green\", "
			"\"-\" lw 1 lc rgb \"blue\"\n"
	      );
	long double m = 0;
	for (int l = 0; l < 3; l++) {
		for (int i = 0; i < n[l]; i++) {
			if (h[l][i][1] > m) m = h[l][i][1];
			printf("\t%Lg\t%Lg\n", h[l][i][0], h[l][i][1]);
		}
		printf("end\n");
	}
}

//void dump_histogram(long double (*h)[2], int n)
//{
//	long double (*a)[2] = xmalloc(n*sizeof*a);
//	memcpy(a, h, n*sizeof*a);
//	accumulate_histogram(a, n);
//	printf("set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
//	printf("set yrange [0:]\n");
//	printf("set format y \"\"\n");
//	printf("unset key\n");
//	printf("plot \"-\" w impulses, \"-\" w lines\n");
//	long double m[4] = {0};
//	for (int i = 0; i < n; i++)
//	{
//		if (h[i][1] > m) m = h[i][1];
//		printf("\t%Lg\t%Lg\n", h[i][0], h[i][1]);
//	}
//	printf("end\n");
//	long double f = m/a[n-1][1];
//	//fprintf(stderr, "m = %Lg\n", m);
//	//fprintf(stderr, "f = %Lg\n", f);
//	//fprintf(stderr, "a = %Lg\n", a[n-1][1]);
//	//fprintf(stderr, "n = %d\n", n);
//	for (int i = 0; i < n; i++)
//		printf("\t%Lg\t%Lg\n", a[i][0], f*a[i][1]);
//	printf("end\n");
//	free(a);
//}

#define RRR 0.333
#define GGG 0.333
#define BBB 0.333

int main(int c, char *v[])
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		//                         0   1
		return EXIT_FAILURE;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	int w, h, pd, nh[4];
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *tmp = xmalloc(w*h*sizeof*tmp);
	if (pd != 3) fail("requires a rgb image");
	long double (*his[4])[2];
	for (int l = 0; l < 4; l++)
		his[l] = xmalloc(w*h*sizeof*his[l]);
	for (int l = 0; l < 3; l++) {
		for (int i = 0; i < w*h; i++)
			tmp[i] = x[3*i+l];
		nh[l] = fill_histogram(his[l], tmp, w*h);
	}
	for (int i = 0; i < w*h; i++)
		tmp[i] = RRR*x[3*i]+GGG*x[3*i+1]+BBB*x[3*i+2];
	nh[3] = fill_histogram(his[3], tmp, w*h);

	dump_histograms(his, nh);

	return EXIT_SUCCESS;
}
