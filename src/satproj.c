// simulator of affine satellite acquisitions
// INPUT:
// 	1. height map
// 	2. texture map
// 	3. projection matrix
// 	4. desired output size
//
// OUTPUT:
// 	1. rendered texture

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bicubic.c"

static void invert_projection(double iA[8], double A[8])
{
	// 
	//  / i \     / a  b \     / l \     / r \         / p \
	// |     | = |        | * |     | + |     | * h + |     |
	//  \ j /     \ c  d /     \ t /     \ s /         \ q /
	//
	double a = A[0], b = A[1], r = A[2], p = A[3];
	double c = A[4], d = A[5], s = A[6], q = A[7];

	double det = a * d - b * c;
	double ia =  d / det;
	double ib = -b / det;
	double ic = -c / det;
	double id =  a / det;
	double ir = -(ia * r + ib * s);
	double is = -(ic * r + id * s);
	double ip = -(ia * p + ib * q);
	double iq = -(ic * p + id * q);

	iA[0] = ia; iA[1] = ib; iA[2] = ir; iA[3] = ip;
	iA[4] = ic; iA[5] = id; iA[6] = is; iA[7] = iq;
}

static void apply_projection(double y[3], double A[8], double x[3])
{
	y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2] + A[3];
	y[1] = A[4] * x[0] + A[5] * x[1] + A[6] * x[2] + A[7];
	y[2] = x[2];
}

void satproj(float *out, int ow, int oh,
		double P[8], float *heights, float *colors,
		int w, int h, int pd)
{
	// P = projection
	// L = localisation
	//
	double L[6];
	invert_projection(L, P);

	fprintf(stderr, "satproj %d %d => %d %d (%d)\n", w, h, ow, oh, pd);
	fprintf(stderr, "projection P = %g %g %g %g  %g %g %g %g\n",
			P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
	fprintf(stderr, "projection L = %g %g %g %g  %g %g %g %g\n",
			L[0], L[1], L[2], L[3], L[4], L[5], L[6], L[7]);

	// fill each point of the output image with the appropriate color
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		double ij[3] = {i, j, 0}, rs[3] = {i, j, 0};
		apply_projection(rs, L, ij);
		float *to = out + pd * (ow * j + i);
		bicubic_interpolation(to, colors, w, h, pd, rs[0], rs[1]);
	}
}






#define MAIN_SATPROJ
#ifdef MAIN_SATPROJ
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
#include "parsenumbers.c"


static void read_n_doubles_from_string(double *out, char *string, int n)
{
	for (int i = 0; i < n; i++)
		out[i] = 0;

	int no;
	double *buf = NULL;
	FILE *f = fopen(string, "r");
	if (f) {
		buf = read_ascii_doubles(f, &no);
		fclose(f);
	} else {
		buf = alloc_parse_doubles(n, string, &no);
	}

	if (no > n) no = n;
	for (int i = 0; i < no; i++)
		out[i] = buf[i];
	free(buf);
}

int main(int c, char *v[])
{
	// input arguments
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s heights.tiff colors.png P.txt ow oh out.png\n", *v);
		//        0 1            2          3     4  5  6
		return 1;
	}
	char *fname_heights = v[1];
	char *fname_colors  = v[2];
	char *fname_pmatrix = v[3];
	int out_w = atoi(v[4]);
	int out_h = atoi(v[5]);
	char *fname_output  = v[6];

	// read input images and matrices
	int w[2], h[2], pd;
	float *heights = iio_read_image_float(fname_heights, w, h);
	float *colors  = iio_read_image_float_vec(fname_colors, w+1, h+1, &pd);
	double P[8] = {0};
	read_n_doubles_from_string(P, fname_pmatrix, 8);

	fprintf(stderr, "heights %d %d\n", w[0], h[0]);
	fprintf(stderr, "colors %d %d %d \n", w[1], h[1], pd);

	// allocate space for output
	float *out = xmalloc(out_w * out_h * pd * sizeof*out);

	// run simulator
	satproj(out, out_w, out_h, P, heights, colors, w[1], h[1], pd);

	// save output
	iio_save_image_float_vec(fname_output, out, out_w, out_h, pd);

	// cleanup and exit
	free(heights);
	free(colors);
	free(out);
	return 0;
}
#endif//MAIN_SATPROJ
