#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fail.c"
#include "xmalloc.c"
#include "random.c"

#include "smapa.h"

static int good_modulus(int n, int p)
{
	if (p < 1)
		fail("bad modulus %d\n", p);

	int r;
	if (n >= 0)
		r = n % p;
	else {
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	if (r < 0 || r >= p)
		fprintf(stderr, "mod(%d,%d)=%d\n",n,p,r);
	assert(r >= 0);
	assert(r < p);
	return r;
}

static void addsample_periodic(float *x, int w, int h, int i, int j, float inc)
{
	int gi = good_modulus(i, w);
	int gj = good_modulus(j, h);
	x[gj*w+gi] += inc;
}

SMART_PARAMETER_SILENT(SEED,0)

static void build_k(float *x, int w, int np, float var, int nc)
{
	srand(SEED());
	float p[np][2], c[2] = {0, 0};
	for (int i = 0; i < w*w; i++)
		x[i] = 0;
	for (int i = 0; i < np; i++)
	{
		p[i][0] = var*random_normal();
		p[i][1] = var*random_normal();
		c[0] += p[i][0] / np;
		c[1] += p[i][1] / np;
	}
	for (int i = 0; i < np; i++)
		addsample_periodic(x, w, w, p[i][0]-c[0], p[i][1]-c[1], 1);
	addsample_periodic(x, w, w, 0, 0, nc);
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 6 && c != 5) {
		fprintf(stderr, "usage:\n\t%s kw np var nc [out]\n", *v);
		//                         0  1  2  3   4   5
		return EXIT_FAILURE;
	}
	int kw = atoi(v[1]);
	int np = atoi(v[2]);
	float var = atof(v[3]);
	int nc = atoi(v[4]);
	char *filename_out = c > 5 ? v[5] : "-";

	int w = 2*kw + 1;
	float *x = xmalloc(w * w * sizeof*x);
	build_k(x, w, np, var, nc);
	iio_write_image_float(filename_out, x, w, w);
	return EXIT_SUCCESS;
}
