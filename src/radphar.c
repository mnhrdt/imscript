// radial p-harmonic function

#include <math.h>

static void radhar(float *X, int w, int h,
		float r1, float r2, float c1, float c2,
		float x0, float y0)
{
	float A = (c2 - c1)/log(r1/r2);
	float B = c1 + A*log(r1);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	{
		float x = i - x0;
		float y = j - y0;
		float r = hypot(x, y);
		float u = -log(r);
		X[j*w+i] = A*u + B;
	}
}

static void radphar(float *X, int w, int h, float p,
		float r1, float r2, float c1, float c2,
		float x0, float y0)
{
	if (p == 2) {
		radhar(X, w, h, r1, r2, c1, c2, x0, y0);
		return;
	}
	float alpha = isfinite(p) ? (p - 2)/(p - 1) : 1;
	float A = (c2 - c1)/(pow(r2,alpha) - pow(r1,alpha));
	float B = c1 - A*pow(r1,alpha);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	{
		float x = i - x0;
		float y = j - y0;
		float r = hypot(x, y);
		float u = pow(r, alpha);
		X[j*w+i] = A*u + B;
	}
}

#define MAIN_RADPHAR

#ifdef MAIN_RADPHAR
#include <stdio.h>
#include "iio.h"
#include "xmalloc.c"
int main(int c, char *v[])
{
	if (c != 10 && c != 11) {
		fprintf(stderr, "usage:\n\t"
		"%s w h p r1 r2 c1 c2 x0 y0 [out.png]\n", *v);
	//        0 1 2 3 4  5  6  7  8  9   10
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	float p = atof(v[3]);
	float r1 = atof(v[4]);
	float r2 = atof(v[5]);
	float c1 = atof(v[6]);
	float c2 = atof(v[7]);
	float x0 = atof(v[8]);
	float y0 = atof(v[9]);
	char *filename_out = c > 10 ? v[10] : "-";

	float *x = xmalloc(w*h*sizeof*x);
	radphar(x, w, h, p, r1, r2, c1, c2, x0, y0);
	iio_save_image_float(filename_out, x, w, h);
	free(x);
	return 0;
}
#endif//MAIN_RADPHAR
