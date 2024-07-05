#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "iio.h"

#include "xmalloc.c"
#include "smapa.h"

SMART_PARAMETER_SILENT(HISMOOTH,0)
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
	if (HISMOOTH() > 0)
		smooth_histogram_rw(h, r, HISMOOTH());
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
	//printf("set format y \"\"\n");
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

static void dump_histogram_c(long double (*h)[2], int n, char *c)
{
	long double (*a)[2] = xmalloc(n*sizeof*a);
	memcpy(a, h, n*sizeof*a);
	accumulate_histogram(a, n);
	printf("set xrange [%Lg:%Lg]\n", h[0][0], h[n-1][0]);
	printf("set yrange [0:]\n");
	printf("unset key\n");
	printf("%s\n", c);
	printf("plot \"-\" w impulses title \"histogram\", \"-\" w lines title \"accumulated histogram\"");
	printf("\n");
	long double m = 0;
	for (int i = 0; i < n; i++)
	{
		if (h[i][1] > m) m = h[i][1];
		printf("\t%Lg\t%Lg\n", h[i][0], h[i][1]);
	}
	printf("end\n");
	long double f = m/a[n-1][1];
	for (int i = 0; i < n; i++)
		printf("\t%Lg\t%Lg\n", a[i][0], f*a[i][1]);
	printf("end\n");
	free(a);
}

static char *help_string_name     = "ghisto";
static char *help_string_version  = "ghisto 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "compute the histogram of an image, in gnuplot format";
static char *help_string_usage    = "usage:\n\t"
"ghisto [-p] [img.png] > histo.g";
static char *help_string_long     =
"Ghisto computes the histogram of an image.\n"
"\n"
"The histogram is printed in a format that can be piped directly to gnuplot.\n"
"Notice that no quantization is made by this program, if the input image\n"
"is floating point, it is likely that all pixel values will be different\n"
"and the histogram will look flat.\n"
"\n"
"Usage: ghisto img.png > histo.g\n"
"   or: cat img.png | ghisto > histo.g\n"
"\n"
"Options:\n"
" -p\t\twrite a png-producing gnuplot program\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Environment:\n"
" HISMOOTH    filter the histogram by a rectangular kernel of this width\n"
" SHOWSTATS   show various statistics (curves, etc)\n"
"\n"
"Examples:\n"
" ghisto img.png | gnuplot                View histogram in a gnuplot window.\n"
" ghisto -p img.png | gnuplot > hist.png  Create png image of the histogram\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
#include "pickopt.c"
int main_ghisto(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);

	bool term_png = pick_option(&c, &v, "p", NULL);
	char *clean_o = pick_option(&c, &v, "c", "");
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		//                         0   1
		return EXIT_FAILURE;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	if (term_png) printf("set term pngcairo\n");

	long double (*his)[2] = xmalloc(w * h * sizeof*his);
	int nh = fill_histogram(his, x, w*h);
	if (*clean_o)
		dump_histogram_c(his, nh, clean_o);
	else
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
