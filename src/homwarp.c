
// y = H(x)
static void apply_homography(double y[2], double H[9], double x[2])
{
	double X = H[0]*x[0] + H[1]*x[1] + H[2];
	double Y = H[3]*x[0] + H[4]*x[1] + H[5];
	double Z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = X / Z;
	y[1] = Y / Z;
}

// compute the inverse homography (inverse of a 3x3 matrix)
static double invert_homography(double invH[9], double H[9])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double det = H[0]*H[4]*H[8] + H[2]*H[3]*H[7] + H[1]*H[5]*H[6]
		   - H[2]*H[4]*H[6] - H[1]*H[3]*H[8] - H[0]*H[5]*H[7];
	invH[0] = ( H[4] * H[8] - H[5] * H[7] ) / det;
	invH[1] = ( H[2] * H[7] - H[1] * H[8] ) / det;
	invH[2] = ( H[1] * H[5] - H[2] * H[4] ) / det;
	invH[3] = ( H[5] * H[6] - H[3] * H[8] ) / det;
	invH[4] = ( H[0] * H[8] - H[2] * H[6] ) / det;
	invH[5] = ( H[2] * H[3] - H[0] * H[5] ) / det;
	invH[6] = ( H[3] * H[7] - H[4] * H[6] ) / det;
	invH[7] = ( H[1] * H[6] - H[0] * H[7] ) / det;
	invH[8] = ( H[0] * H[4] - H[1] * H[3] ) / det;
	return det;
}


#include <math.h>
#include "bicubic_gray.c"
#include "bilinear_interpolation.c"

static
float nearest_neighbor_interpolator(float *x, int w, int h, float p, float q)
{
	int ip = round(p);
	int iq = round(q);
	if (ip < 0) ip = 0;
	if (iq < 0) iq = 0;
	if (ip >= w) ip = w - 1;
	if (iq >= h) iq = h - 1;
	return x[w*iq+ip];
}

typedef float (*gray_interpolator_t)(float *,int,int,float,float);

void homwarp(float *X, int W, int H, double M[9], float *x, int w, int h, int o)
{
	gray_interpolator_t u = bicubic_interpolation_gray;
	if (o == 0) u = nearest_neighbor_interpolator;
	if (o == 2) u = bilinear_interpolation_at;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		double p[2] = {i, j}, Mp[2];
		apply_homography(Mp, M, p);
		X[j*W+i] = u(x, w, h, Mp[0], Mp[1]);
	}
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
#include "parsenumbers.c"
#include "pickopt.c"
int main(int c, char *v[])
{
	int order = atoi(pick_option(&c, &v, "i", "3"));
	if (c < 4 || c > 6)
		return fprintf(stderr, "usage:\n\t"
				"%s [-i {0|2|3}] hom w h [in [out]]\n", *v);
		//                0                1   2 3  4   5
	double H[9];
	read_n_doubles_from_string(H, v[1], 9);
	int ow = atoi(v[2]);
	int oh = atoi(v[3]);
	char *filename_in  = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = xmalloc(ow * oh * pd * sizeof*y);

	for (int i = 0; i < pd; i++)
		homwarp(y + i*ow*oh, ow, oh, H, x + i*w*h, w, h, order);

	iio_save_image_float_split(filename_out, y, ow, oh, pd);
	return 0;
}
