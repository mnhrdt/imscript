#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define xmalloc malloc
#define FOR(i,n) for(int i=0;i<(n);i++)

int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	float (*y)[2*w][pd] = xmalloc(4*w*h*pd*sizeof(float));

	FOR(j,h) FOR(i,w) FOR(l,pd)
	{
		float g = x[(j*w+i)*pd + l];
		int si = 2*w-i-1;
		int sj = 2*h-j-1;
		y [ j] [ i] [l] = g;
		y [sj] [ i] [l] = g;
		y [ j] [si] [l] = g;
		y [sj] [si] [l] = g;
	}

	iio_write_image_float_vec(out, y[0][0], 2*w, 2*h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
