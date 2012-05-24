#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"

#include "xmalloc.c"

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
	return r;
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
	printf("set yrange [0:]\n");
	printf("set format y \"\"\n");
	printf("unset key\n");
	printf("plot \"-\" w impulses, \"-\" w lines\n");
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

int main(int c, char *v[])
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
