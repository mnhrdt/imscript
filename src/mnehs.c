#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "iio.h"

#include "xmalloc.c"
#include "getpixel.c"
#include "bicubic.c"

#define DONT_USE_TEST_MAIN
#include "rpc.c"

static void apply_projection(double y[3], double P[8], double x[3])
{
	y[0] = P[0] * x[0] + P[1] * x[1] + P[2] * x[2] + P[3];
	y[1] = P[4] * x[0] + P[5] * x[1] + P[6] * x[2] + P[7];
	y[2] = x[2];
}


float eval_using_k33(float *x, int w, int h, int i, int j, float k[9])
{
	getsample_operator p = getsample_2;
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
	float f = 0.125;
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

void mnehs_rpc_warp(float *warp,
		float *h, int wh, int hh, int pd,
		float *a, int wa, int ha,
		float *b, int wb, int hb,
		struct rpc *ra, struct rpc *rb)
{
	for (int j = 0; j < hh; j++)
	for (int i = 0; i < wh; i++)
	{
		float vh = getsample_nan(h, wh, hh, 1, i, j, 0);
		double ijh[3] = {i, j, vh};
		double paijh[3], pbijh[3];

		eval_rpci(paijh, ra, i, j, vh);
		eval_rpci(pbijh, rb, i, j, vh);
		fprintf(stderr, "(%d %d %g) => (%g %g)\n", i, j, vh, paijh[0], paijh[1]);

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
SMART_PARAMETER(NWARPS,1)

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

static float laplacian8_at(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_0;

	float r = -8 * p(x, w, h, i  , j  )
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1)
		     + p(x, w, h, i+1, j+1)
		     + p(x, w, h, i+1, j-1)
		     + p(x, w, h, i-1, j-1)
		     + p(x, w, h, i-1, j+1);

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
		bicubic_interpolation_nans(&va, a, wa, ha, 1, paijh[0], paijh[1]);
		bicubic_interpolation_nans(&vb, b, wb, hb, 1, pbijh[0], pbijh[1]);
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
			if (!isfinite(amb[ij])) continue;

			float ax = laplacian_at(h, ow, oh, i, j);
			ax -= Q[ij] * (Q[ij] * h[ij] + amb[ij]) / alpha2;
			h[ij] += TAU() * ax;
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


static void center_projection(double P[8], double cx, double cy)
{
	// point "c" must be fixed at h = 0
	P[3] = cx - P[0]*cx - P[1]*cy;
	P[7] = cy - P[4]*cx - P[5]*cy;
}

#include <stdbool.h>
#include "pickopt.c"


int main_warp(int c, char *v[])
{
	// input arguments
	bool do_center = pick_option(&c, &v, "c", NULL);
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png Pa.txt Pb.txt in.tiff amb.png\n", *v);
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

	fprintf(stderr, "pa = %g %g %g %g %g %g %g %g\n",
		PA[0], PA[1], PA[2], PA[3], PA[4], PA[5], PA[6], PA[7]);
	fprintf(stderr, "pb = %g %g %g %g %g %g %g %g\n",
		PB[0], PB[1], PB[2], PB[3], PB[4], PB[5], PB[6], PB[7]);

	// perform centering, if necessary
	if (do_center) {
		fprintf(stderr, "centering projection matrices\n");
		center_projection(PA, wa/2, ha/2);
		center_projection(PB, wb/2, hb/2);
	}

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

int main_warprpc(int c, char *v[])
{
	// input arguments
	bool do_center = pick_option(&c, &v, "c", NULL);
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png Pa.rpc Pb.rpc in.tiff amb.png\n", *v);
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

	//double PA[8], PB[8];
	//read_n_doubles_from_string(PA, matrix_pa, 8);
	//read_n_doubles_from_string(PB, matrix_pb, 8);
	struct rpc rpca[1];
	struct rpc rpcb[1];
	read_rpc_file_xml(rpca, matrix_pa);
	read_rpc_file_xml(rpcb, matrix_pb);

	if (pd != pdb)
		fail("input pair has different color depth");


	// allocate space for output image
	float *out = xmalloc(wi * hi * pd * sizeof*out);

	// run the algorithm
	mnehs_rpc_warp(out, h0, wi,hi,pd, a,wa,ha, b,wb,hb, rpca, rpcb);

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
	bool do_only_warp = pick_option(&c, &v, "w", NULL);
	if (do_only_warp) return main_warp(c, v);
	bool do_center = pick_option(&c, &v, "c", NULL);
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

	// perform centering, if necessary
	if (do_center) {
		center_projection(PA, wa/2, ha/2);
		center_projection(PB, wb/2, hb/2);
	}

	// allocate space for output image
	float *out = xmalloc(wi * hi * sizeof*out);

	// run the algorithm
	float alpha2 = ALPHA()*ALPHA();
	int niter = NITER();
	int nwarps = NWARPS();
	for (int i = 0; i < nwarps; i++)
	{
		mnehs_affine(out, h0, wi,hi, a,wa,ha, b,wb,hb, PA, PB, alpha2, niter);
		memcpy(h0, out, wi * hi * sizeof*h0);
	}

	// save the output image
	iio_save_image_float(filename_out, out, wi, hi);

	// cleanup and exit
	free(out);
	free(h0);
	free(a);
	free(b);
	return 0;
}
int main(int c,char*v[]){return main_warp(c,v);}
#endif//MAIN_MNEHS
