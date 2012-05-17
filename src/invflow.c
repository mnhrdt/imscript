#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "bicubic.c"



static void flowinv_iter(float *v, float *u, int w, int h)
{
	float (*V)[w][2] = (void*)v;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p[2];
		float qx = i + V[j][i][0];
		float qy = j + V[j][i][1];
		bicubic_interpolation(p, u, w, h, 2, qx, qy);
		V[j][i][0] = -p[0];
		V[j][i][1] = -p[1];
	}
}

static void flowinv(float *v, float *u, int w, int h, int niter, int epsil)
{
	fprintf(stderr, "flowinv %d %d\n", w, h);
	//float *v = vo;
	//float *t = xmalloc(2*w*h*sizeof*t);

	for (int i = 0; i < w*h*2; i++)
		v[i] = -u[i];

	for (int i = 0; i < niter; i++)
	//{
		flowinv_iter(v, u, w, h);
		//void *tmp = t; t = v; v = tmp;
	//}

	//if (v != vo)
	//	for (int i = 0; i < w*h*2; i++)
	//		vo[i] = t[i];

	//free(t);
}

int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s niter epsil [in [out]]\n", *v);
		//                          0 1     2      3   4
		return EXIT_FAILURE;
	}
	int niter = atoi(v[1]);
	float epsil = atof(v[2]);
	char *infile = c > 3 ? v[3] : "-";
	char *outfile = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) fail("2D vector field expected");
	float *y = xmalloc(2*w*h*sizeof*y);
	flowinv(y, x, w, h, niter, epsil);
	iio_save_image_float_vec(outfile, y, w, h, 2);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
