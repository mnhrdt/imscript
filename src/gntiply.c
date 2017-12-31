#include <assert.h>
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

static void pix_weighted_acc(float *y, float *x, float w, int p)
{
	for (int i = 0; i < p; i++)
		y[i] += w * x[i];
}

static void pix_weighted_sum(float *y, float **x, float *w, int p, int n)
{
	for (int k = 0; k < p; k++)
		y[k] = 0;
	for (int k = 0; k < p; k++)
	for (int i = 0; i < n; i++)
		y[k] += w[i] * x[i][k];
}

void ntiply_epx2(float *y, float *x, int w, int h, int p)
{
	// Eric's Pixel eXpansion algorithm, by Eric Johnston of LucasArts
	//
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
		float *t[4] = {
			pix_get(y, 2*w, 2*h, p, 2*i  , 2*j  ),
			pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j  ),
			pix_get(y, 2*w, 2*h, p, 2*i  , 2*j+1),
			pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j+1),
		};
		for (int k = 0; k < 4; k++)
			pix_cp(t[k], P, p);
#define eq(a,b) pix_eq(a,b,p)
		if (eq(C,A) && !eq(C,D) && !eq(A,B)) pix_cp(t[0], A, p);
		if (eq(A,B) && !eq(A,C) && !eq(B,D)) pix_cp(t[1], B, p);
		if (eq(D,C) && !eq(D,B) && !eq(C,A)) pix_cp(t[2], C, p);
		if (eq(B,D) && !eq(B,A) && !eq(D,C)) pix_cp(t[3], D, p);
#undef eq
	}
}

void ntiply_epx3(float *y, float *x, int w, int h, int p)
{
// pseudocode from https://en.wikipedia.org/wiki/Pixel-art_scaling_algorithms
//
// A B C --\  1 2 3
// D E F    > 4 5 6
// G H I --/  7 8 9
//
// 1=E; 2=E; 3=E; 4=E; 5=E; 6=E; 7=E; 8=E; 9=E;
// IF D==B AND D!=H AND B!=F :
//         1=D
// IF (D==B AND D!=H AND B!=F AND E!=C) OR (B==F AND B!=D AND F!=H AND E!=A) :
//         2=B
// IF B==F AND B!=D AND F!=H :
//         3=F
// IF (H==D AND H!=F AND D!=B AND E!=A) OR (D==B AND D!=H AND B!=F AND E!=G) :
//         4=D
// 5=E
// IF (B==F AND B!=D AND F!=H AND E!=I) OR (F==H AND F!=B AND H!=D AND E!=C) :
//         6=F
// IF H==D AND H!=F AND D!=B
//         7=D
// IF (F==H AND F!=B AND H!=D AND E!=G) OR (H==D AND H!=F AND D!=B AND E!=I) :
//         8=H
// IF F==H AND F!=B AND H!=D :
//         9=F
//
// Note: the C code below is obtained from this pseudocode by search&replace

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float *A = pix_get(x, w, h, p, i-1, j-1);
		float *B = pix_get(x, w, h, p, i  , j-1);
		float *C = pix_get(x, w, h, p, i+1, j-1);
		float *D = pix_get(x, w, h, p, i-1, j  );
		float *E = pix_get(x, w, h, p, i  , j  );
		float *F = pix_get(x, w, h, p, i+1, j  );
		float *G = pix_get(x, w, h, p, i-1, j+1);
		float *H = pix_get(x, w, h, p, i  , j+1);
		float *I = pix_get(x, w, h, p, i+1, j+1);
		float *t[9] = {
			pix_get(y, 3*w, 3*h, p, 3*i-1, 3*j-1),
			pix_get(y, 3*w, 3*h, p, 3*i  , 3*j-1),
			pix_get(y, 3*w, 3*h, p, 3*i+1, 3*j-1),
			pix_get(y, 3*w, 3*h, p, 3*i-1, 3*j  ),
			pix_get(y, 3*w, 3*h, p, 3*i  , 3*j  ),
			pix_get(y, 3*w, 3*h, p, 3*i+1, 3*j  ),
			pix_get(y, 3*w, 3*h, p, 3*i-1, 3*j+1),
			pix_get(y, 3*w, 3*h, p, 3*i  , 3*j+1),
			pix_get(y, 3*w, 3*h, p, 3*i+1, 3*j+1),
		};
		for (int k = 0; k < 9; k++)
			pix_cp(t[k], E, p);
#define eq(a,b) pix_eq(a,b,p)
		if ( eq(D,B) && !eq(D,H) && !eq(B,F) )
			pix_cp(t[0], D, p);
		if ( (eq(D,B) && !eq(D,H) && !eq(B,F) && !eq(E,C))
			|| (eq(B,F) && !eq(B,D) && !eq(F,H) && !eq(E,A)) )
			pix_cp(t[1], B, p);
		if ( eq(B,F) && !eq(B,D) && !eq(F,H) )
			pix_cp(t[2], F, p);
		if ( (eq(H,D) && !eq(H,F) && !eq(D,B) && !eq(E,A))
			|| (eq(D,B) && !eq(D,H) && !eq(B,F) && !eq(E,G)) )
			pix_cp(t[3], D, p);
		//pix_cp(t[4], E, p);
		if ( (eq(B,F) && !eq(B,D) && !eq(F,H) && !eq(E,I))
			|| (eq(F,H) && !eq(F,B) && !eq(H,D) && !eq(E,C)) )
			pix_cp(t[5], F, p);
		if ( eq(H,D) && !eq(H,F) && !eq(D,B) )
			pix_cp(t[6], D, p);
		if ( (eq(F,H) && !eq(F,B) && !eq(H,D) && !eq(E,G))
			|| (eq(H,D) && !eq(H,F) && !eq(D,B) && !eq(E,I)) )
			pix_cp(t[7], H, p);
		if (eq(F,H) && !eq(F,B) && !eq(H,D))
			pix_cp(t[8], F, p);
#undef eq
	}
}

static void hq2x_fill_lut(float t[4*256])
{
	for (int i = 0; i < 4*256; i++)
		t[i] = 0.25;
}

void ntiply_hq2x(float *y, float *x, int w, int h, int p)
{
	// build lookup table (could be cached if necessary)
	float t[4*256];
	hq2x_fill_lut(t);

	// initialize large image to zero (to be accumulated into)
	for (int i = 0; i < 4*w*h*p; i++)
		y[i] = 0;

	// traverse the small image
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		// extract pixel values at 3x3 neighborhood
		float *A[9] = {
			pix_get(x, w, h, p, i-1, j-1), // 0
			pix_get(x, w, h, p, i  , j-1), // 1
			pix_get(x, w, h, p, i+1, j-1), // 2
			pix_get(x, w, h, p, i-1, j  ), // 3
			pix_get(x, w, h, p, i  , j  ), // 4 <- center
			pix_get(x, w, h, p, i+1, j  ), // 5
			pix_get(x, w, h, p, i-1, j+1), // 6
			pix_get(x, w, h, p, i  , j+1), // 7
			pix_get(x, w, h, p, i+1, j+1), // 8
		};

		// local binary pattern (8 bits)
		int lbp = 0;
		for (int b = 0; b < 8; b++)
			lbp = 2 * lbp + pix_eq(A[b], A[4], p);
		assert(lbp >= 0 && lbp < 256);

		// 0 1 2
		// 3 4 5
		// 6 7 8
		float *p00[4] = {A[0], A[1], A[3], A[4]};
		float *p10[4] = {A[1], A[2], A[4], A[5]};
		float *p01[4] = {A[3], A[4], A[6], A[7]};
		float *p11[4] = {A[4], A[5], A[7], A[8]};
		float *q00 = pix_get(y, 2*w, 2*h, p, 2*i  , 2*j  );
		float *q10 = pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j  );
		float *q01 = pix_get(y, 2*w, 2*h, p, 2*i  , 2*j+1);
		float *q11 = pix_get(y, 2*w, 2*h, p, 2*i+1, 2*j+1);
		pix_weighted_acc(q00, *p00, 4*t[4*lbp + 0], p);
		pix_weighted_acc(q10, *p10, 3*t[4*lbp + 1], p);
		pix_weighted_acc(q01, *p01, 2*t[4*lbp + 2], p);
		pix_weighted_acc(q11, *p11, 1*t[4*lbp + 3], p);
	}
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
	else if (0 == strcmp(m, "epx") && f == 3)
		ntiply_epx3(y, x, w, h, p);
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
