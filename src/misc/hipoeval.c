
void hypoeval(float *o, float *v, float *s, int w, int h, int n)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int p = s[j*w+i];
		if (p < 0) p = 0;
		if (p >= n) p = n - 1;
		o[j*w+i] = v[(j*w+i)*n+p];
	}
}

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 4) return fprintf(stderr, "usage:\n\t"
			"%s vol slice out\n", *v);
	//                0 1   2     3
	char *filename_vol = v[1];
	char *filename_sli = v[2];
	char *filename_out = v[3];
	int W, H, pd, w, h;
	float *vol = iio_read_image_float_vec(filename_vol, &W, &H, &pd);
	assert(pd == 256);
	float *sli = iio_read_image_float(filename_sli, &w, &h);
	assert(w == W && h == H);
	float *out = malloc(w*h*sizeof*out);
	hypoeval(out, vol, sli, w, h, pd);
	iio_write_image_float_vec(filename_out, out, w, h, 1);
	return 0;
}
