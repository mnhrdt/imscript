#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "xmalloc.c"

static int bound(int x, int min, int max)
{
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

static void crop(float *out, int *cw, int *ch, float *in, int w, int h, int pd,
		int x0, int y0, int xf, int yf)
{
	if (xf < 0) xf = w - xf;
	if (yf < 0) yf = h - yf;
	x0 = bound(x0, 0, w);
	xf = bound(xf, 0, w);
	y0 = bound(y0, 0, h);
	yf = bound(yf, 0, h);
	if (x0 >= xf) fail("bad crop x");
	if (y0 >= yf) fail("bad crop y");

	*cw = xf - x0;
	*ch = yf - y0;
	float (*x)[w][pd] = (void*)in;
	float (*y)[*cw][pd] = (void*)out;

	for (int j = 0; j < *ch; j++)
	for (int i = 0; i < *cw; i++)
	for (int l = 0; l < pd; l++)
		y[j][i][l] = x[j+y0][i+x0][l];
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c < 5 || c > 7) {
		fprintf(stderr, "usage:\n\t%s x0 y0 xf yf [in [out]]\n", *v);
		//                          0 1  2  3  4   5   6
		return EXIT_FAILURE;
	}
	int x0 = atoi(v[1]);
	int y0 = atoi(v[2]);
	int xf = atoi(v[3]);
	int yf = atoi(v[4]);
	char *filename_in = c > 5 ? v[5] : "-";
	char *filename_out = c > 6 ? v[6] : "-";

	int w, h, pd;
	float *image_in = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *image_out = xmalloc(w*h*pd*sizeof*image_out);

	int cw, ch;
	crop(image_out, &cw, &ch, image_in, w, h, pd,
			x0, y0, xf, yf);

	iio_save_image_float_vec(filename_out, image_out, cw, ch, pd);
	return EXIT_SUCCESS;
}
