// invert a vector field

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "bicubic.c"

static bool checkbounds(int a, int x, int b)
{
	return a <= x && x < b;
}

static void flowinv_init(float *v, float *u, int w, int h)
{
	float (*U)[w][2] = (void*)u;
	float (*V)[w][2] = (void*)v;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < 2; l++)
		V[j][i][l] = -U[j][i][l];

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float px = i + U[j][i][0];
		float py = j + U[j][i][1];
		int ipx = px;
		int ipy = py;
		if (checkbounds(0, ipx, w) && checkbounds(0, ipy, h))
			for (int l = 0; l < 2; l++)
				V[ipy][ipx][l] = -U[j][i][l];
		if (checkbounds(0, ipx+1, w) && checkbounds(0, ipy, h))
			for (int l = 0; l < 2; l++)
				V[ipy][ipx+1][l] = -U[j][i][l];
		if (checkbounds(0, ipx, w) && checkbounds(0, ipy+1, h))
			for (int l = 0; l < 2; l++)
				V[ipy+1][ipx][l] = -U[j][i][l];
		if (checkbounds(0, ipx+1, w) && checkbounds(0, ipy+1, h))
			for (int l = 0; l < 2; l++)
				V[ipy+1][ipx+1][l] = -U[j][i][l];
	}
}

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

static float flowinv_eval(float *v, float *u, int w, int h)
{
	static int idcount = 0;
	float (*U)[w][2] = (void*)u;
	float (*V)[w][2] = (void*)v;

	float (*sav)[w][2] = xmalloc(2 * w * h * sizeof(float));
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < 2; l++)
		sav[j][i][l] = 0;

	int ppt = (w+h)/30;

	long double r = 0;
	int cx = 0;
	for (int j = ppt; j < h-ppt; j++)
	for (int i = ppt; i < w-ppt; i++)
	{
		float p[2];
		float qx = i + V[j][i][0];
		float qy = j + V[j][i][1];
		bicubic_interpolation(p, u, w, h, 2, qx, qy);
		r += hypot(V[j][i][0] + p[0], V[j][i][1] + p[1]);
		sav[j][i][0] = V[j][i][0] + p[0];
		sav[j][i][1] = V[j][i][1] + p[1];
		cx += 1;
	}
	idcount += 1;
	char fname[0x100];
	snprintf(fname, 0x100, "/tmp/fresi_%03d_%s.flo",
			idcount/2, idcount%2?"odd":"even");
	iio_save_image_float_vec(fname, sav[0][0], w, h, 2);
	free(sav);
	return r/cx;
}

static void flowinv_printeval(float *v, float *u, int w, int h)
{
	float r1 = flowinv_eval(v, u, w, h);
	float r2 = flowinv_eval(u, v, w, h);
	fprintf(stderr, "e %g\t%g\n", r1, r2);
}

static void flowinv(float *v, float *u, int w, int h, int niter, int epsil)
{
	for (int i = 0; i < w*h*2; i++)
		v[i] = -u[i];
	if (epsil > 0)
		flowinv_init(v, u, w, h);

	flowinv_printeval(v, u, w, h);
	for (int i = 0; i < niter; i++)
	{
		flowinv_iter(v, u, w, h);
		flowinv_printeval(v, u, w, h);
	}
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
