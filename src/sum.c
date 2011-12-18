#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define BAD_MAX(a,b) (a)<(b)?(b):(a)


static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "xmalloc: out of memory\n"));
	return new;
}

typedef float (*extension_operator_float)(float*,int,int,int,int,int,int);

static float extend_float_image_by_zero(float *xx, int w, int h, int pd,
		int i, int j, int l)
{
	float (*x)[w][pd] = (void*)xx;
	if (i < 0 || j < 0 || i > w-1 || j > h-1 || l < 0 || l > pd-1)
		return 0;
	else
		return x[j][i][l];
}

int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s in1 [in2 [out]]\n", *v);
		//                          0 1    2    3
		return EXIT_FAILURE;
	}
	char *infile1 = v[1];
	char *infile2 = c > 2 ? v[2] : "-";
	char *outfile = c > 3 ? v[3] : "-";

	int w[3], h[3], pd[3];
	float *x1 = iio_read_image_float_vec(infile1, w+0, h+0, pd+0);
	float *x2 = iio_read_image_float_vec(infile2, w+1, h+1, pd+1);

	w[2] = BAD_MAX(w[0], w[1]);
	h[2] = BAD_MAX(h[0], h[1]);
	pd[2] = BAD_MAX(pd[0], pd[1]);

	float (*y)[w[2]][pd[2]] = xmalloc(w[2] * h[2] * pd[2] * sizeof(float));

	for (int i = 0; i < w[2] * h[2] * pd[2]; i++)
		((float*)(y))[i]= 0;

	extension_operator_float p = extend_float_image_by_zero;
	for (int j = 0; j < h[2]; j++)
	for (int i = 0; i < w[2]; i++)
	for (int l = 0; l < pd[2]; l++)
		y[j][i][l] = p(x1,w[0],h[0],pd[0], i,j,l)
				+ p(x2,w[1],h[1],pd[1], i,j,l);

	iio_save_image_float_vec(outfile, y[0][0], w[2], h[2], pd[2]);
	return EXIT_SUCCESS;
}
