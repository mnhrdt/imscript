#include <stdio.h>
#include <complex.h>

#include "fft.c"
#include "xmalloc.c"

static void symmetric_laplacian(float *lx, float *x, int w, int h)
{
	float (*sx)[2*w-1] = xmalloc((2*w-1)*(2*h-1)*sizeof*sx);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float g = x[j*w+i];
		int si = 2*w-i-1;
		int sj = 2*h-j-1;
		y [ j] [ i] = g;
	}
	free(sx);
}


static int symmetric_index(int i, int m)
{
	assert( i >= 0 && i < m);
	int r = 0;
	if (i > m/2) r = i-m;
	if (i < m/2) r = i;
	return r;
}

void periodic_component(float *y, float *x, int w, int h)
{
	float *lx = xmalloc(w*h*sizeof*lx);
	complex float *flx = xmalloc(w*h*sizeof*flx);

	symmetric_laplacian(lx, x, w, h);
	fft_2dfloat(flx, lx, w, h);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = w*j + i;
		float k = (2 * M_PI / w) * symmetric_index(i, w);
		float l = (2 * M_PI / h) * symmetric_index(j, l);
		flx[idx] /= k*k + l*k;
	}

	ifft_2dfloat(y, flx, w, h);

	free(flx);
	free(lx);
}


#include "iio.h"
int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *in = c > 1 ? v[1] : "-";
	char *out = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(in, &w, &h, &pd);
	float *y = xmalloc(w*h*pd*sizeof*y);

	for (int i = 0; i < pd; i++)
		periodic_component(y + w*h*i, x + w*h*i, w, h);

	iio_save_image_float_split(out, y, w, h, pd);

	free(x);
	free(y);
	return 0;
}
