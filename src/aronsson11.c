// Aronsson infinity-harmonic function

#include <math.h>

static void aronsson11(float *X, int w, int h, float alpha, float x0, float y0)
{
	alpha *= 4*atan(1)/180;
	float c = cos(alpha);
	float s = sin(alpha);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	{
		float x = i - x0;
		float y = j - y0;
		float rotx = c*x - s*y;
		float roty = s*x + c*y;
		float u = fabs(rotx) - fabs(roty);
		X[j*w+i] = u;
	}
}

#define MAIN_ARONSSON11

#ifdef MAIN_ARONSSON11
#include <stdio.h>
#include "iio.h"
#include "xmalloc.c"
int main(int c, char *v[])
{
	if (c != 6 && c != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s w h alpha0 x0 y0 [out.png]\n", *v);
	//        0 1 2 3      4  5   6
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	float alpha = atof(v[3]);
	float x0 = atof(v[4]);
	float y0 = atof(v[5]);
	char *filename_out = c > 6 ? v[6] : "-";

	float *x = xmalloc(w*h*sizeof*x);
	aronsson11(x, w, h, alpha, x0, y0);
	iio_save_image_float(filename_out, x, w, h);
	free(x);
	return 0;
}
#endif//MAIN_ARONSSON11
