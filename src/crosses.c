// make a transparent image with crosses at indicated positions


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define CROSS_SIDE 50
#define CROSS_HOLE 23

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

void fill_cross_back(uint8_t *xx, int w, int h, int i, int j, bool up)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return;
	uint8_t (*x)[w][4] = (void*)xx;
	int k = up ? 0 : 0xff;
	x[j][i][0] = k;
	x[j][i][1] = k;
	x[j][i][2] = k;
	x[j][i][3] = 0x80;
}

void fill_cross_rgba(uint8_t *x, int w, int h, int a, int b)
{
	if (a < 0 || b < 0 || a >= w || b >= h)
		return;
	//for (int i = a-CROSS_SIDE; i <= a+CROSS_SIDE; i++)
	//	fill_cross_point(x, w, h, i, b);
	//for (int i = b-CROSS_SIDE; i <= b+CROSS_SIDE; i++)
	//	fill_cross_point(x, w, h, a, i);
	for (int j = -CROSS_SIDE; j <= CROSS_SIDE; j++)
	for (int i = -CROSS_SIDE; i <= CROSS_SIDE; i++)
	{
		float n = hypot(i, j);
		if (!j || !i)
			fill_cross_point(x, w, h, a+i, b+j);
		else if (n < CROSS_SIDE && n > CROSS_HOLE)
			fill_cross_back(x, w, h, a+i, b+j, (i+j)*(i-j) > 0);
	}
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
	if (c < 4 && c%2) {
		fprintf(stderr, "usage:\n\t%s out w h [p0x p0y [p1x ... ]", *v);
		//                          0 1   2 3  4   5    6
			return EXIT_FAILURE;
	}
	for (int i = 0 ; i < c; i++)
		fprintf(stderr, "crosses ARGV[%d] = %s\n", i, v[i]);
	char *filename_out = v[1];;
	int w = atoi(v[2]);
	int h = atoi(v[3]);
	int np = (c-4)/2;
	int p[2*np];
	for (int i = 0; i < 2*np; i++)
		p[i] = atoi(v[4+i]);

	uint8_t *x = xmalloc(w * h * 4);
	fill_crosses_rgba(x, w, h, p, np);

	iio_write_image_uint8_vec(filename_out, x, w, h, 4);
	free(x);

	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
