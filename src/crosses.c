// make a transparent image with crosses at indicated positions


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define CROSS_SIDE 20

void fill_cross_point(uint8_t *xx, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return;
	uint8_t (*x)[w][4] = (void*)xx;
	int k = (i+j)%2 ? 0 : 0xff;
	x[j][i][0] = k;
	x[j][i][1] = k;
	x[j][i][2] = k;
	x[j][i][3] = 0xff;
}

void fill_cross_rgba(uint8_t *x, int w, int h, int a, int b)
{
	for (int i = a-CROSS_SIDE; i <= a+CROSS_SIDE; i++)
		fill_cross_point(x, w, h, i, b);
	for (int i = b-CROSS_SIDE; i <= b+CROSS_SIDE; i++)
		fill_cross_point(x, w, h, a, i);
}

void fill_crosses_rgba(uint8_t *x, int w, int h, int *p, int np)
{
	for (int i = 0; i < w * h * 4; i++)
		x[i] = 0;

	for (int i = 0; i < np; i++)
		fill_cross_rgba(x, w, h, p[2*i], p[2*i+1]);
}


#ifndef OMIT_MAIN

#include "fail.c"
#include "xmalloc.c"
#include "iio.h"

int main(int c, char *v[])
{
	if (c < 3 && 0==c%2) {
		fprintf(stderr, "usage:\n\t%s w h [p0x p0y [p1x p1y ... ]", *v);
		//                          0 1 2  3   4    5   6
			return EXIT_FAILURE;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int np = (c-3)/2;
	int p[2*np];
	for (int i = 0; i < 2*np; i++)
		p[i] = atoi(v[3+i]);

	uint8_t *x = xmalloc(w * h * 4);
	fill_crosses_rgba(x, w, h, p, np);

	iio_save_image_uint8_vec("-", x, w, h, 4);
	free(x);

	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
