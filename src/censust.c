// census transform of a single image

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


static void pack_bits_into_bytes(unsigned char *out, int *bits, int nbits)
{
	int nbytes = ceil(nbits/8.0);
	for (int i = 0; i < nbytes; i++)
	{
		out[i] = 0;
		for (int j = 0; j < 8; j++)
			out[i] = out[i] * 2 + bits[8*i+j];
	}
}

static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	return x[(i+j*w)*pd + l];
}


static void census_at(unsigned char *out, float *x, int w, int h, int pd,
		int winradius, int p, int q)
{
	int side = 2*winradius + 1;
	int nbits = pd * (side * side - 1);
	int bits[nbits];
	int cx = 0;
	for (int l = 0; l < pd; l++)
	for (int j = -winradius; j <= winradius; j++)
	for (int i = -winradius; i <= winradius; i++)
	{
		float a = getsample_nan(x, w, h, pd, p    , q    , l);
		float b = getsample_nan(x, w, h, pd, p + i, q + j, l);
		if (i || j)
			bits[cx++] = a < b;
	}
	assert(cx == nbits);

	pack_bits_into_bytes(out, bits, nbits);
}

static void color_census_transform(unsigned char *y, int opd,
		float *x, int w, int h, int pd, int winradius)
{

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		census_at(y + opd * (w * j + i), x, w, h, pd, winradius, i, j);
}

#include "iio.h"
#include "pickopt.c"
int main_censust(int c, char *v[])
{
	// process command line arguments
	char *radius_opt = pick_option(&c, &v, "r", "1");
	if ((c == 2 && 0 == strcmp("-h", v[1])) || (c != 1 && c != 2 && c != 3))
	{
		fprintf(stderr, "usage:\n\t%s [-r radius] [in [out]]\n", *v);
		//                          0              1   2
		return 1;
	}
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";
	int winradius = atoi(radius_opt);

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	// size of neighborhood in bits and in bytes
	int side = 2 * winradius + 1;
	int nbits = pd * (side * side - 1);
	int nbytes = ceil(nbits / 8.0);
	fprintf(stderr, "output nbits=%d nbytes=%d\n", nbits, nbytes);

	// allocate output image
	unsigned char *y = malloc(w * h * nbytes);

	// compute stuff
	color_census_transform(y, nbytes, x, w, h, pd, winradius);

	// save, cleanup, and exit
	iio_write_image_uint8_vec(filename_out, y, w, h, nbytes);
	free(x);
	free(y);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_censust(c, v); }
#endif
