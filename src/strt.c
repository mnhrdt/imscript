// structure tensor computations



#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "xmalloc.c"


typedef float (*extension_operator_float)(float*, int, int, int, int);

static float extend_float_image_by_zero(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i > w-1 || j > h-1)
		return 0;
	else
		return x[j*w+i];
}

static float extend_float_image_constant(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static double sqr(double x) { return x; }

static void compute_structure_tensor_here(float tensor[3],
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, int w, int h,
		int i, int j)
{
	extension_operator_float p = extend_float_image_constant;

	int n = kside * kside;
	tensor[0] = tensor[1] = tensor[2] = 0;

	for (int k = 0; k < n; k++)
	{
		int ii = i + wo[k][0];
		int jj = j + wo[k][1];
		float wvk = wv ? wv[k] : 1;
		tensor[0] += wvk * sqr(p(gx, w, h, ii, jj));
		tensor[1] += wvk * p(gx,w,h, ii, jj) * p(gy,w,h, ii, jj);
		tensor[2] += wvk * sqr(p(gy, w, h, ii, jj));
	}
}

static void fill_window_offsets(int (*wo)[2], int kside)
{
	int idx = 0;
	int kradius = (kside - 1)/2;
	for (int j = 0; j < kside; j++)
	for (int i = 0; i < kside; i++)
	{
		wo[idx][0] = i - kradius;
		wo[idx][1] = j - kradius;
		idx += 1;
	}
}

static void fill_window_values(float *wv, int (*wo)[2], int kside, float sigma)
{
	int n = kside*kside;
	for (int i = 0; i < n; i++)
	{
		if (sigma < 0)
			wv[i] = 1;
		else {
			float r = hypot(wo[i][0], wo[i][1]);
			wv[i] = exp(-sqr(r/sigma));
		}
	}
	float m = 0;
	for (int i = 0; i < n; i++)
		m += wv[i];
	for (int i = 0; i < n; i++)
		wv[i] /= m;

	//if (0) {
	//	FILE *f = fopen("/tmp/lk.kkk", "w");
	//	fprintf(f, "P2\n%d %d\n65535\n", kside, kside);
	//	for (int i = 0; i < n; i++)
	//		fprintf(f, "%g\n", wv[i]);
	//	fclose(f);
	//}
}

static void setpixel_float_image_vec(float *x, int w, int h, int pd,
		int i, int j, int l, float v)
{
	if (i < 0 || j < 0 || l < 0 || i>=w || j>=h || l >= pd)
		exit(fprintf(stderr, "bad setpixel (%d %d %d)[%d %d %d]\n",
					w, h, pd, i, j, l));
	float (*xx)[w][pd] = (void*)x;
	xx[j][i][l] = v;
}

static void compute_structure_tensor_field_fancy(float *out_st,
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, int w, int h)
{

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float atwa[3];
		compute_structure_tensor_here(atwa, wv, wo, kside, gx, gy, w, h,
				i, j);
		for (int l = 0; l < 3; l++)
			setpixel_float_image_vec(out_st, w,h,3, i,j,l, atwa[l]);
	}
}

// out_st is filled to be an image of 3 channels (the components of the tensor)
static void compute_structure_tensor_field(float *out_st,
		int kside, double sigma,
		float *gx, float *gy, int w, int h)
{
	int wo[kside*kside][2]; fill_window_offsets(wo, kside);
	float wv[kside*kside]; fill_window_values(wv, wo, kside, sigma);
	compute_structure_tensor_field_fancy(out_st, wv, wo, kside, gx,gy, w,h);
}

// "out_stf" is filled to be an image of 7 channels, containing the structure
// tensor T at each point, together with several related descriptors
//
// 0: T11
// 1: T12
// 2: T22
// 3: large eigenvalue
// 4: small eigenvalue
// 5: principal direction x
// 6: principal direction y
//
// Note that the principal direction is determined up to a rotation of 180
// degrees
//
static void compute_structure_tensor_field_ultra_fancy(float *out_stf,
		int x, int w, int h, int kside, double sigma)
{
	int wo[kside*kside][2]; fill_window_offsets(wo, kside);
	float wv[kside*kside]; fill_window_values(wv, wo, kside, sigma);
	float *tt = malloc(3 * w * h * sizeof(float));
	float *gx = malloc(2 * w * h * sizeof(float));
	float *gy = gx + w*h;
	for (int i = 0; i < w*h; i++)
		gx[i] = gy[i] = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int ij = j*w + i;
		int Ij = j*w + i + 1;
		int iJ = j*w + i + w;
		gx[ij] = x[Ij] - x[ij];
		gy[ij] = x[iJ] - x[ij];
	}
	compute_structure_tensor_field_fancy(tt, wv, wo, kside, gx,gy, w,h);
	for (int i = 0; i < w * h; i++)
	{
		double a = tt[3*i+0];
		double b = tt[3*i+1];
		double c = tt[3*i+2];
		double T = a + c;
		double D = a*c - b*b;
		assert(D >= 0);
		assert(T*T - 4*D >= 0);
		double lambda = ( T + sqrt( T*T - 4*D ) ) / 2;
		double mu     = ( T - sqrt( T*T - 4*D ) ) / 2;
		assert(lambda >= mu);
		assert(mu >= 0);
		double vx = 1, vy = 0; // default eigenvector
		if (fabs(b) > fabs(c)) {
			vx = b;
			vy = lambda - a;
		} else {
			vx = lambda - d;
			vy = c;
		}
		out_stf[7*i+0] = a;
		out_stf[7*i+1] = b;
		out_stf[7*i+2] = c;
		out_stf[7*i+3] = lambda;
		out_stf[7*i+4] = mu;
		out_stf[7*i+5] = vx;
		out_stf[7*i+6] = vy;

	}
	free(gx);
	free(tt);
}

#ifndef OMIT_MAIN_STRT
#define MAIN_STRT
#endif//OMIT_MAIN_STRT

#ifdef MAIN_STRT
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s kside ksigma [in [out]]\n", *v);
		//                          0 1     2       3   4
		return 1;
	}
	int kside = atoi(v[1]);
	double ksigma = atof(v[2]);
	char *filename_in  = c >= 3 ? v[3] : "-";
	char *filename_out = c >= 4 ? v[4] : "-";

	int w, h;
	float *x  = iio_read_image_float(filename_in, &w, &h);
	float *t  = xmalloc(w * h * 7 * sizeof*t);

	compute_structure_tensor_field_ultra_fancy(t, x, w, h, kside, ksigma);

	iio_save_image_float_vec(filename_out, t, w, h, 7);

	free(x); free(t);
	return 0;
}
#endif//MAIN_STRT
