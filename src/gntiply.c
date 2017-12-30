#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ntiply_naive(float *y, float *x, int w, int h, int p, int f)
{
	for (int j = 0; j < f*h; j++)
	for (int i = 0; i < f*w; i++)
	for (int l = 0; l < p  ; l++)
		y[p*(f*w*j+i)+l] = x[p*(w*(j/f)+i/f)+l];
}

// getpixel with nearest-neighbour extrapolation
static float *pix_get(float *x, int w, int h, int p, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x + p*(j*w+i);
}

// function to compare whether two pixels are close enough
static bool pix_eq(float *a, float *b, int p)
{
	//double r = 0;
	//for (int i = 0; i < p; i++)
	//	r = hypot(r, a[i] - b[i]);
	//return r < 40;
	for (int i = 0; i < p; i++)
		if (a[i] != b[i])
			return false;
	return true;
}

static void pix_cp(float *y, float *x, int p)
{
	for (int i = 0; i < p; i++)
		y[i] = x[i];
}

void ntiply_epx2(float *y, float *x, int w, int h, int p)
{
	//    A    --\ 1 2
	//  C P B  --/ 3 4
	//    D
	//
	//   1=P; 2=P; 3=P; 4=P;
	//   IF C==A AND C!=D AND A!=B => 1=A
	//   IF A==B AND A!=C AND B!=D => 2=B
	//   IF D==C AND D!=B AND C!=A => 3=C
	//   IF B==D AND B!=A AND D!=C => 4=D

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float *P = pix_get(x, w, h, p, i  , j  );
		float *A = pix_get(x, w, h, p, i  , j-1);
		float *B = pix_get(x, w, h, p, i+1, j  );
		float *C = pix_get(x, w, h, p, i-1, j  );
		float *D = pix_get(x, w, h, p, i  , j+1);
		float *t[4] = { pix_get(y, 2*w, 2*h, p, 2*i  , 2*j  ),
		                pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j  ),
		                pix_get(y, 2*w, 2*h, p, 2*i  , 2*j+1),
		                pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j+1)
		};
		for (int k = 0; k < 4; k++)
			pix_cp(t[k], P, p);
		if (pix_eq(C,A,p) && !pix_eq(C,D,p) && !pix_eq(A,B,p))
			pix_cp(t[0], A, p);
		if (pix_eq(A,B,p) && !pix_eq(A,C,p) && !pix_eq(B,D,p))
			pix_cp(t[1], B, p);
		if (pix_eq(D,C,p) && !pix_eq(D,B,p) && !pix_eq(C,A,p))
			pix_cp(t[2], C, p);
		if (pix_eq(B,D,p) && !pix_eq(B,A,p) && !pix_eq(D,C,p))
			pix_cp(t[3], D, p);
	}
}

void ntiply_hq2x(float *y, float *x, int w, int h, int p)
{
	fprintf(stderr, "WARNING: hq2x not implemented\n");
	ntiply_naive(y, x, w, h, p, 2);
}

void ntiply_generic(
		float *y,  // output image, of size (w*f) * (h*f) * p
		float *x,  // input image, of size w * h * p
		int w,     // width of the input image
		int h,     // height of the input image
		int p,     // pixel dimension
		int f,     // scaling factor
		char *m    // name of the scaling method
		)
{
	if      (0 == strcmp(m, "epx") && f == 2)
		ntiply_epx2(y, x, w, h, p);
	else if (0 == strcmp(m, "hqx") && f == 2)
		ntiply_hq2x(y, x, w, h, p);
	else
		ntiply_naive(y, x, w, h, p, f);
}

#include "iio.h"
int main_gntiply(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3)
		return fprintf(stderr, "usage:\n\t%s m f [in [out]]\n", *v);
		//                                 0 1 2  3   4
	char *method_id    = v[1];
	int scale_factor   = atof(v[2]);
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int f = scale_factor;
	int w, h, p;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &p);
	float *y = malloc(f*f*w*h*p*sizeof*y);

	ntiply_generic(y, x, w, h, p, f, method_id);

	iio_write_image_float_vec(filename_out, y, f*w, f*h, p);
	free(x);
	free(y);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_gntiply(c, v); }
#endif
