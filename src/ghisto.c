#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "iio.h"

#include "xmalloc.c"
#include "smapa.h"

SMART_PARAMETER_SILENT(HISMOTTH,0)
SMART_PARAMETER_SILENT(SHOWSTATS,0)

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

static int fill_histogram(long double (*h)[2], float *in_x, int n)
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

static void accumulate_histogram(long double (*h)[2], int n)
{
	for (int i = 1; i < n; i++)
		h[i][1] += h[i-1][1];
}

static void print_gnuplot_stats_string(long double (*)[2], int, float);

static void dump_histogram(long double (*h)[2], int n)
{
	long double (*a)[2] = xmalloc(n*sizeof*a);
	memcpy(a, h, n*sizeof*a);
	accumulate_histogram(a, n);
	printf("set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
	printf("set yrange [0:]\n");
	printf("set format y \"\"\n");
	if (SHOWSTATS() > 0) {
		printf("set samples 1000\nset key left\n");
		printf("N(m,s,x)=exp(-(x-m)**2/(2*s*s))/(s*sqrt(2*pi))\n");
		printf("L(m,s,x)=exp(-sqrt(2)*abs(x-m)/s)/(s*sqrt(2))\n");
		printf("U0(a,b,x)=x<a?0:(x>b?0:1/(b-a))\n");
		printf("U(m,s,x)=U0(m-s*sqrt(3),m+s*sqrt(3),x)\n");
		//printf("C(m,s,x)=1/((1+((x-m)/s)**2)*(pi*s))\n");
	}
	else
		printf("unset key\n");
	printf("plot \"-\" w impulses title \"histogram\", \"-\" w lines title \"accumulated histogram\"");
	if (SHOWSTATS() > 0) {
		print_gnuplot_stats_string(h, n, 1/SHOWSTATS());
	}
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

int main_ghisto(int c, char *v[])
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		//                         0   1
		return EXIT_FAILURE;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	long double (*his)[2] = xmalloc(w * h * sizeof*his);
	int nh = fill_histogram(his, x, w*h);
	dump_histogram(his, nh);

	free(x);
	free(his);
	return EXIT_SUCCESS;
}

static void print_gnuplot_stats_string(long double (*h)[2], int n, float bs)
{
	long double mass = 0, avg = 0, var = 0;
	for (int i = 0; i < n; i++)
		if (isfinite(h[i][1]))
			mass += h[i][1];
	for (int i = 0; i < n; i++)
		if (isfinite(h[i][1]))
			avg += h[i][0] * h[i][1];
	avg /= mass;
	for (int i = 0; i < n; i++)
		if (isfinite(h[i][1]))
			var += (h[i][0] - avg) * (h[i][0] - avg) * h[i][1];
	var /= mass;
	fprintf(stderr, "mass = %g\n", (double)mass);
	fprintf(stderr, "nmass = %g\n", (double)mass/n);
	fprintf(stderr, "mu = %g\n", (double)avg);
	fprintf(stderr, "sigma = %g\n", sqrt(var));
	printf(", %g*N(%g,%g,x) w lines title \"gaussian\"",
			(double)mass*bs, (double)avg, sqrt(var));
	printf(", %g*L(%g,%g,x) w lines title \"laplacian\"",
			(double)mass*bs, (double)avg, sqrt(var));
	printf(", %g*U(%g,%g,x) w lines title \"uniform\"",
			(double)mass*bs, (double)avg, sqrt(var));
	//printf(", %g*C(%g,%g,x) w lines title \"cauchy\"",
	//		(double)mass*bs, (double)avg, sqrt(var));
	//note: Cauchy estimation is much hairier ("avg" and "var" do not work)
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_ghisto(c, v); }
#endif
