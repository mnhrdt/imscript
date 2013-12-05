// from a set of pairs to an "inpainting" mask

#include <math.h>

void pairsinp(float *o, int w, int h, float *p, int np)
{
	for (int i = 0; i < 2*w*h; i++)
		o[i] = NAN;
	for (int k = 0; k < np; k++)
	{
		float from[2], to[2], d[2];
		int ifrom[2];
		for (int l = 0; l < 2; l++)
		{
			from[l] = p[4*k+l];
			to[l] = p[4*k+2+l];
			ifrom[l] = round(from[l]);
			d[l] = to[l] - from[l];
		}
		int i = ifrom[0];
		int j = ifrom[1];
		if (i >= 0 && i < w && j >= 0 && j < h)
		{
			int idx = w*j + i;
			for (int l = 0; l < 2; l++)
				o[2*idx+l] = d[l];
		}
	}
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "xmalloc.c"
#include "parsenumbers.c"

int main(int c, char **v)
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s w h <pairs.txt >img.tiff\n", *v);
		return 0;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);

	int n;
	float *p = read_ascii_floats(stdin, &n);
	n /= 4;

	float *x = xmalloc(2*w*h*sizeof*x);
	pairsinp(x, w, h, p, n);

	iio_save_image_float_vec("-", x, w, h, 2);

	free(x);
	free(p);
	return 0;
}
