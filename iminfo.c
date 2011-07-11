#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"

#include "fragments.c"
#include "statistics.c"

static void printvals(float *x, int n)
{
	for (int i = 0; i < n; i++)
	{
		printf("%g", x[i]);
		if (i < n - 1)
			putchar(' ');
	}
}


int main(int c, char *v[])
{
	if (c != 1 && c != 2) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		return EXIT_FAILURE;
	}
	char *filename = c > 1 ? v[1] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename, &w, &h, &pd);
	printf("image containing %dx%d %d-dimensional pixels\n", w, h, pd);
	printf("first pixel: "); printvals(x, pd); putchar('\n');
	printf("central pixel: "); printvals(x+pd*((w*h)/2), pd); putchar('\n');
	struct statistics_float s;
	statistics_getf(&s, x, w*h*pd);
	print_stats(stdout, &s, "samples");
	if (pd > 1) for (int i = 0; i < pd; i++) {
		char buf[100]; snprintf(buf, 100, "%dth component", i);
		statistics_getf_stride(&s, x + i, w*h, pd);
		print_stats(stdout, &s, buf);
	}
	return EXIT_SUCCESS;
}
