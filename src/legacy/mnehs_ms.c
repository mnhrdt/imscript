#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

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

#include "smapa.h"

SMART_PARAMETER(TAU,0.25)
SMART_PARAMETER(ALPHA,1)
SMART_PARAMETER(NITER,1)
SMART_PARAMETER(NWARPS,1)
SMART_PARAMETER(NSCALES,1)

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

static void nusavec(char *fmt, int idx, float *x, int w, int h, int pd)
{
	char fname[FILENAME_MAX];
	snprintf(fname, FILENAME_MAX, fmt, idx);
	iio_write_image_float_vec(fname, x, w, h, pd);
}

static int global_scale;

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

	nusavec("Q_%d.tiff", global_scale, Q, ow, oh, 1);
	nusavec("amb_%d.tiff", global_scale, amb, ow, oh, 1);

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
	nusavec("h_%d.tiff", global_scale, h, ow, oh, 1);

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

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
//	if (PREFILTER() > 0) {
//		float *fin = xmalloc(iw*ih*sizeof*fin);
//		for (int j = 0; j < ih; j++)
//		for (int i = 0; i < iw; i++)
//		{
//			//if (isfinite(in[iw*j+i])) {
//			//	float a[9], m = 0;
//			//	for (int ii=-1;ii<=1;ii++)
//			//	for (int jj=-1;jj<=1;jj++)
//			//		a[3*ii+jj] = p(in, iw, ih, i+ii, j+jj);
//			//	int cx = 0;
//			//	for (int k = 0; k < 9; k++)
//			//		if (isfinite(a[k])) {
//			//			m += a[k];
//			//			cx += 1;
//			//		}
//			//	fin[iw*j + i] = cx ? m/cx : NAN;
//			//} else {
//				float a[4], m = 0;
//				a[0] = p(in, iw, ih, i+1, j);
//				a[1] = p(in, iw, ih, i-1, j);
//				a[2] = p(in, iw, ih, i, j+1);
//				a[3] = p(in, iw, ih, i, j-1);
//				int cx = 0;
//				for (int k = 0; k < 4; k++)
//					if (isfinite(a[k])) {
//						m += a[k];
//						cx += 1;
//					}
//				fin[iw*j + i] = cx ? m/cx : NAN;
//			//}
//		}
//		in = fin;
//	}
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i, 2*j);
		a[1] = p(in, iw, ih, 2*i+1, 2*j);
		a[2] = p(in, iw, ih, 2*i, 2*j+1);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
//	if (PREFILTER() > 0) free(in);
}

// evaluate a bilinear cell at the given point
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = p(in, iw, ih, round((i-0.5)/2), round((j-0.5)/2));
		//out[ow*j+i] = bilinear_interpolation(in, iw, ih, (i-0.5)/2, (j-0.5)/2);
}

void mnehs_affine_ms(float *out, float *in, int w, int h,
		float *a, int aw, int ah,
		float *b, int bw, int bh,
		double PA[8], double PB[8],
		float alpha2, int niter, int scale)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1)
	{
		int  ws = ceil( w/2.0); int  hs = ceil( h/2.0);
		int aws = ceil(aw/2.0); int ahs = ceil(ah/2.0);
		int bws = ceil(bw/2.0); int bhs = ceil(bh/2.0);
		float *is = xmalloc( ws *  hs * sizeof*is);
		float *os = xmalloc( ws *  hs * sizeof*os);
		float *as = xmalloc(aws * bhs * sizeof*as);
		float *bs = xmalloc(bws * bhs * sizeof*bs);

		double PAs[8] = { PA[0], PA[1], PA[2]/2, PA[3]/2,
			          PA[4], PA[5], PA[6]/2, PA[7]/2 };
		double PBs[8] = { PB[0], PB[1], PB[2]/2, PB[3]/2,
			          PB[4], PB[5], PB[6]/2, PB[7]/2 };
		zoom_out_by_factor_two(is,  ws,  hs, in, w,  h);
		zoom_out_by_factor_two(as, aws, ahs, a, aw, ah);
		zoom_out_by_factor_two(bs, bws, bhs, b, bw, bh);
		mnehs_affine_ms(os, is, ws, hs, as, aws, ahs, bs, bws, bhs,
				PAs, PBs, alpha2, niter, scale - 1);
		zoom_in_by_factor_two(init, w, h, os, ws, hs);

		free(is);
		free(os);
		free(as);
		free(bs);
	} else {
		for (int i = 0; i < w * h; i++)
		       init[i] = in[i];	
	}
	global_scale = scale;
	mnehs_affine(out, init, w,h, a,aw,ah, b,bw,bh, PA, PB, alpha2, niter);
	free(init);
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
	iio_write_image_float_vec(filename_out, out, wi, hi, pd);

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
	int nscales = NSCALES();
	mnehs_affine_ms(out, h0, wi, hi, a, wa, ha, b, wb, hb, PA, PB, alpha2, niter, nscales);
	//for (int i = 0; i < nwarps; i++)
	//{
	//	mnehs_affine(out, h0, wi,hi, a,wa,ha, b,wb,hb, PA, PB, alpha2, niter);
	//	memcpy(h0, out, wi * hi * sizeof*h0);
	//}

	// save the output image
	iio_write_image_float(filename_out, out, wi, hi);

	// cleanup and exit
	free(out);
	free(h0);
	free(a);
	free(b);
	return 0;
}
int main(int c,char*v[]){return main_compute(c,v);}
#endif//MAIN_MNEHS
