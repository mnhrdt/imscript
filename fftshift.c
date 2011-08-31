#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"


#define xmalloc malloc
#define FOR(i,n) for(int i=0;i<(n);i++)


static int good_modulus(int n, int p)
{
	assert(p);
	if (p < 0) return good_modulus(n, -p);

	int r;
	if (n >= 0)
		r = n % p;
	else
	{
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	assert(r >= 0);
	assert(r < p);
	return r;
}

typedef float (*extension_operator_float)(float*,int,int,int,int,int,int);

static float extend_float_image_periodic(float *xx, int w, int h, int pd,
		int i, int j, int l)
{
	float (*x)[w][pd] = (void*)xx;
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	l = good_modulus(l, pd);
	return x[j][i][l];
}

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
	float (*y)[w][pd] = xmalloc(w*h*pd*sizeof(float));

	extension_operator_float p = extend_float_image_periodic;

	FOR(j,h) FOR(i,w) FOR(l,pd)
		y[j][i][l] = p(x, w,h,pd, i+w/2, j+h/2, l);

	iio_save_image_float_vec(out, y[0][0], w, h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
