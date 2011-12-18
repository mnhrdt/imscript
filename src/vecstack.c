#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define xmalloc malloc

int main(int c, char *v[])
{
	if (c < 3) {
		fprintf(stderr, "usage:\n\t%s [chan1 ...] > out\n", *v);
		return EXIT_FAILURE;
	}
	float *chan[c-1];
	int w[c-1], h[c-1];
	for (int i = 0; i < c-1; i++)
		chan[i] = iio_read_image_float(v[i+1], w + i, h + i);
	float (*out)[c-1] = xmalloc(*w * *h * (c-1) * sizeof*out);
	for (int i = 0; i < c-1; i++) {
		if (w[i] != *w || h[i] != *w)
			fprintf(stderr,"warning: %dth image size mismatch\n",i);
		else
			for (int j = 0; j < *w * *h; j++)
				out[j][i] = chan[i][j];
	}
	iio_save_image_float_vec("-", out[0], *w, *h, c-1);
	return EXIT_SUCCESS;
}
