#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"

#include "xmalloc.c"
#include "getpixel.c"
#include "bicubic.c"

static void apply_projection(double y[3], double P[8], double x[3])
{
	y[0] = P[0] * x[0] + P[1] * x[1] + P[2] * x[2] + P[3];
	y[1] = P[4] * x[0] + P[5] * x[1] + P[6] * x[2] + P[7];
	y[2] = x[2];
}


float eval_using_k33(float *x, int w, int h, int i, int j, float k[9])
{
	getsample_operator p = getsample_1;
	int cx = 0;
	float r = 0;
	for (int jj = j-1; jj <= j+1; jj++)
	for (int ii = i-1; ii <= i+1; ii++)
		r += k[cx++] * p(x, w, h, 1, ii, jj, 0);
	return r;
}

static void fill_gradient(float *g, float *x, int w, int h)
{
	float kdx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
	float kdy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
	float f = 0.25;
	//float kdx[9] = {0, 0, 0,  0, -1, 1,  0, 0, 0};
	//float kdy[9] = {0, 0, 0,  0, -1, 0,  0, 1, 0};
	//float f = 1;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		g[2*(j*w+i)+0]  = f*eval_using_k33(x, w, h, i, j, kdx);
		g[2*(j*w+i)+1]  = f*eval_using_k33(x, w, h, i, j, kdy);
	}
}

#include "bicubic.c"


void mnehs_affine_warp(float *warp,
		float *h, int wh, int hh, int pd,
		float *a, int wa, int ha,
		float *b, int wb, int hb,
		double PA[8], double PB[8])
{
	for (int j = 0; j < hh; j++)
	for (int i = 0; i < wh; i++)
	{
		float vh = getsample_nan(h, wh, hh, 1, i, j, 0);
		double ijh[3] = {i, j, vh};
		double paijh[3], pbijh[3];
		apply_projection(paijh, PA, ijh);
		apply_projection(pbijh, PB, ijh);
		float va[pd], vb[pd];
		bicubic_interpolation(va, a, wa, ha, pd, paijh[0], paijh[1]);
		bicubic_interpolation(vb, b, wb, hb, pd, pbijh[0], pbijh[1]);
		for (int l = 0; l < pd; l++)
			warp[pd*(j*wh+i)+l] = va[l] - vb[l];
	}
}

#include "smapa.h"

SMART_PARAMETER(TAU,0.25)
SMART_PARAMETER(ALPHA,1)
SMART_PARAMETER(NITER,1)

static float laplacian_at(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = -4 * p(x, w, h, i  , j  )
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1);

	return r;
}

static float hbar_at(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = 0
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1);

	return r/4;
}

// Modèle Numérique d'Élévation Horn Schunck (cas affine)
void mnehs_affine(float *out_h, float *init_h, int ow, int oh,
		float *a, int wa, int ha,
		float *b, int wb, int hb,
		double PA[8], double PB[8],
		float alpha2, int niter)
{
	// allocate temporary images
	float *h   = xmalloc(ow * oh * sizeof*h);     // h-increment
	float *Q   = xmalloc(ow * oh * sizeof*Q);      // Q
	float *amb = xmalloc(ow * oh * sizeof*amb);    // A-B (warped)
	float *ga  = xmalloc(2 * wa * ha * sizeof*ga); // grad(A)
	float *gb  = xmalloc(2 * wb * hb * sizeof*gb); // grad(B)

	// gradient of A and B
	fill_gradient(ga, a, wa, ha);
	fill_gradient(gb, b, wb, hb);

	// fill images q and amb (a minus b)
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float h0 = getsample_nan(init_h, ow, oh, 1, i, j, 0);
		double ijh[3] = {i, j, h0}, paijh[3], pbijh[3];
		apply_projection(paijh, PA, ijh);
		apply_projection(pbijh, PB, ijh);
		float va, vb, vga[2], vgb[2];
		bicubic_interpolation(&va, a, wa, ha, 1, paijh[0], paijh[1]);
		bicubic_interpolation(&vb, b, wb, hb, 1, pbijh[0], pbijh[1]);
		bicubic_interpolation(vga, ga, wa, ha, 2, paijh[0], paijh[1]);
		bicubic_interpolation(vgb, gb, wb, hb, 2, pbijh[0], pbijh[1]);
		float gapa = vga[0] * PA[2] + vga[1] * PA[6];
		float gbpb = vgb[0] * PB[2] + vgb[1] * PB[6];
		Q  [j*ow+i] = gapa - gbpb;
		amb[j*ow+i] = va - vb;
	}

	iio_save_image_float_vec("Q.tiff", Q, ow, oh, 1);
	iio_save_image_float_vec("amb.tiff", amb, ow, oh, 1);

	// initialize h
	for (int i = 0; i < ow * oh; i++)
		h[i] = 0;

	// run the iterations (without warps)
	for (int iter = 0; iter < niter; iter++)
	{
		for (int j = 0; j < oh; j++)
		for (int i = 0; i < ow; i++)
		{
			int ij = j * ow + i;

			//float hbar = hbar_at(h, ow, oh, i, j);
			//float ax = (-hbar * Q[ij] - amb[ij]) * Q[ij];
			//ax /= alpha2 + Q[ij] * Q[ij];
			//h[ij] = hbar + ax;
			//float ax = 0;

			float ax = laplacian_at(h, ow, oh, i, j);
			////ax += laplacian_at(init_h, ow, oh, i, j);
			ax -= Q[ij] * (Q[ij] * h[ij] + amb[ij]) / alpha2;
			////ax -= ( amb[idx]) / alpha2;
			h[ij] += TAU() * tanh(ax);
		}
	}
	iio_save_image_float_vec("h.tiff", h, ow, oh, 1);

	// update result
	for (int i = 0; i < ow * oh; i++)
		out_h[i] = init_h[i] + h[i];
	
	// cleanup and exit
	free(h);
	free(Q);
	free(amb);
	free(ga);
	free(gb);
}

#define MAIN_MNEHS
#ifdef MAIN_MNEHS
#include <stdio.h>
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

int main_warp(int c, char *v[])
{
	// input arguments
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png Pa.txt Pb.txt in.tiff out.tiff\n", *v);
		//        0  1     2     3      4     5       6
		return 1;
	}
	char *filename_a   = v[1];
	char *filename_b   = v[2];
	char *matrix_pa    = v[3];
	char *matrix_pb    = v[4];
	char *filename_in  = v[5];
	char *filename_out = v[6];

	// read input images and matrices
	int wa, wb, wi, ha, hb, hi, pd, pdb;
	float *a  = iio_read_image_float_vec(filename_a, &wa, &ha, &pd);
	float *b  = iio_read_image_float_vec(filename_b, &wb, &hb, &pdb);
	float *h0 = iio_read_image_float(filename_in, &wi, &hi);
	double PA[8], PB[8];
	read_n_doubles_from_string(PA, matrix_pa, 8);
	read_n_doubles_from_string(PB, matrix_pb, 8);
	if (pd != pdb)
		fail("input pair has different color depth");

	// allocate space for output image
	float *out = xmalloc(wi * hi * pd * sizeof*out);

	// run the algorithm
	mnehs_affine_warp(out, h0, wi,hi,pd, a,wa,ha, b,wb,hb, PA, PB);

	// save the output image
	iio_save_image_float_vec(filename_out, out, wi, hi, pd);

	// cleanup and exit
	free(out);
	free(h0);
	free(a);
	free(b);
	return 0;
}

int main_compute(int c, char *v[])
{
	// input arguments
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png Pa.txt Pb.txt in.tiff out.tiff\n", *v);
		//        0  1     2     3      4     5       6
		return 1;
	}
	char *filename_a   = v[1];
	char *filename_b   = v[2];
	char *matrix_pa    = v[3];
	char *matrix_pb    = v[4];
	char *filename_in  = v[5];
	char *filename_out = v[6];

	// read input images and matrices
	int wa, wb, wi, ha, hb, hi;
	float *a  = iio_read_image_float(filename_a, &wa, &ha);
	float *b  = iio_read_image_float(filename_b, &wb, &hb);
	float *h0 = iio_read_image_float(filename_in, &wi, &hi);
	double PA[8], PB[8];
	read_n_doubles_from_string(PA, matrix_pa, 8);
	read_n_doubles_from_string(PB, matrix_pb, 8);

	// allocate space for output image
	float *out = xmalloc(wi * hi * sizeof*out);

	// run the algorithm
	float alpha2 = ALPHA()*ALPHA();
	int niter = NITER();
	mnehs_affine(out, h0, wi,hi, a,wa,ha, b,wb,hb, PA, PB, alpha2, niter);

	// save the output image
	iio_save_image_float(filename_out, out, wi, hi);

	// cleanup and exit
	free(out);
	free(h0);
	free(a);
	free(b);
	return 0;
}
int main(int c,char*v[]){return main_compute(c,v);}
#endif//MAIN_MNEHS
