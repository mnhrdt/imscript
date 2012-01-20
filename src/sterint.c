#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"


#include "xmalloc.c"

static void fill_line(float *y, float *x, int n)
{
	// compute hole locations
	int hole[n][2], nholes = 0;
	bool inhole = false;
	for (int i = 0; i < n; i++) {
		if (inhole && !isnan(x[i])) { // end of hole
			assert(i > 0 && isnan(x[i-1]));
			fprintf(stderr, "\twhile in hole %d, found a number (%g) at %d\n", nholes, x[i], i);
			hole[nholes-1][1] = i-1;
			inhole = false;
		} else if (!inhole && isnan(x[i])) { // start of hole
			if (i > 0) assert(!isnan(x[i-1]));
			fprintf(stderr, "\tout of a hole, found nan at %d, starting h=%d\n", i, nholes);
			hole[nholes][0] = i;
			inhole = true;
			nholes += 1;
		}
		y[i] = x[i];
	}
	if (inhole)
		hole[nholes-1][1] = n-1;

	fprintf(stderr, "got %d holes\n", nholes);

	// fill-in holes
	for (int h = 0; h < nholes; h++)
	{
		int a = hole[h][0];
		int b = hole[h][1];
		fprintf(stderr, "hole %d: [%d %d]\n", h, a, b);
		assert(a <= b);
		if (a > 0) assert(!isnan(x[a-1]) && isnan(x[a]));
		if (b < n-1) assert(!isnan(x[b+1]) && isnan(x[b]));
		float first = NAN, last = NAN;
		if (a > 0) first = x[a-1];
		if (b < n-1) last = x[b+1];
		if (isnan(first) && !isnan(last))
			first = last;
		if (isnan(last) && !isnan(first))
			last = first;
		fprintf(stderr, "f, l = %g, %g\n", first, last);
		if (isnan(first)) assert(isnan(last));
		if (!isnan(first)) {
			float alpha = (last - first)/(2 + b - a);
			float beta = first -  alpha * (a - 1);
			for (int i = a; i <= b; i++)
				y[i] = alpha * i + beta;
		}
	}
}

void stereo_interpolation(float *y, float *x, int w, int h)
{
	for (int i = 0; i < h; i++)
		fill_line(y + i*w, x + i*w, w);

	for (int i = 0; i < w; i++) {
		float cin[h], cout[h];
		for (int j = 0; j < h; j++)
			cin[j] = y[w*j+i];
		fill_line(cout, cin, h);
		for (int j = 0; j < h; j++)
			y[w*j+i] = cout[j];
	}
}

int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h;
	float *x = iio_read_image_float(in, &w, &h);
	float *y = xmalloc(w*h*sizeof*y);
	stereo_interpolation(y, x, w, h);
	iio_save_image_float(out, y, w, h);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
