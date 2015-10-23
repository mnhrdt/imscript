#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include "iio.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define OMIT_MAIN_STRT
#include "strt.c"

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// cut a line with a segment (returns true if they cut)
static bool cut_line_with_segment(double out[2], double line[3],
		double p[2], double q[2])
{
	// points in "oriented" projective coordinates
	double pp[3] = {p[0], p[1], 1};
	double qq[3] = {q[0], q[1], 1};

	// sign of each point (says on which side of the line each point is)
	double sp = scalar_product(pp, line);
	double sq = scalar_product(qq, line);

	// if signs are different, the line crosses the segment
	if (sp * sq < 0) {
		// line trough points p and q
		double pq[3]; vector_product(pq, pp, qq);

		// intersection of "line" and "pq"
		double ii[3]; vector_product(ii, pq, line);

		// recover affine coordinates
		out[0] = ii[0] / ii[2];
		out[1] = ii[1] / ii[2];
		return true;
	}
	return false;
}

static bool cut_line_with_rectangle(double out_a[2], double out_b[2],
		double line[3], double rec_from[2], double rec_to[4])
{
	// four vertices of the rectangle
	double v[4][2] = {
		{ rec_from[0], rec_from[1] },
		{ rec_to[0]  , rec_from[1] },
		{ rec_to[0]  , rec_to[1]   },
		{ rec_from[0], rec_to[1]   }
	};

	// intersections with each of the edges
	bool xP[4]; // whether it intersects
	double x[4][2]; // where it intersects
	for (int i = 0; i < 4; i++)
		xP[i] = cut_line_with_segment(x[i], line, v[i], v[ (i+1)%4 ] );

	// write output
	int n_intersections = xP[0] + xP[1] + xP[2] + xP[3];
	if (n_intersections == 2) { // generic case: 2 intersections
		int cx = 0;
		for (int i = 0; i < 4; i++)
			if (xP[i])
			{
				double *out = cx ? out_b : out_a;
				out[0] = x[i][0];
				out[1] = x[i][1];
				cx += 1;
			}
		return true;
	}
	return false;
}

static void cut_two_lines(double out_p[2], double l[3], double m[3])
{
	// compute crossing point in projective coordinates
	double p[3];
	vector_product(p, l, m);

	// recover affine coordinates
	out_p[0] = p[0] / p[2];
	out_p[1] = p[1] / p[2];
}

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

struct tdip_state;
typedef bool local_orientation_t(double [2], struct tdip_state *, int, int);

// data structure to allow an uniform access to different orientation methods
struct tdip_state {
	// input and output data
	int w, h; // size of the input image
	float *x; // image data
	int tside;
	double aradius;
	float *transform;

	// options
	int nrandom;
	int ntensor;
	local_orientation_t *o;

	// temporary data, used when necessary
	float *structure_tensor;
};

static
bool get_positive_gradient(double g[2], float *x, int w, int h, int i, int j)
{
	if (!insideP(w, h, i, j)) return false;
	if (!insideP(w, h, i+1, j+1)) return false;
	double xij = x[ (j + 0)*w + (i + 0) ]; if(!isfinite(xij))return false;
	double xIj = x[ (j + 0)*w + (i + 1) ]; if(!isfinite(xIj))return false;
	double xiJ = x[ (j + 1)*w + (i + 0) ]; if(!isfinite(xiJ))return false;
	double xIJ = x[ (j + 1)*w + (i + 1) ]; if(!isfinite(xIJ))return false;
	//if (xij <= 9 || xIj <= 9 || xiJ <= 9 || xIJ <= 9) return false;
	g[0] = 0.5 * (xIj - xij + xIJ - xiJ);
	g[1] = 0.5 * (xiJ - xij + xIJ - xIj);
	//g[0] = xIj - xij;
	//g[1] = xiJ - xij;
	return true;
}

static bool local_orientation_by_gradient(double g[2], struct tdip_state *e,
		int i, int j)
{
	int w = e->w;
	int h = e->h;
	float *x = e->x;
	if (!insideP(w, h, i, j)) return false;
	if (!insideP(w, h, i+1, j+1)) return false;
	double xij = x[ (j + 0)*w + (i + 0) ]; if(!isfinite(xij))return false;
	double xIj = x[ (j + 0)*w + (i + 1) ]; if(!isfinite(xIj))return false;
	double xiJ = x[ (j + 1)*w + (i + 0) ]; if(!isfinite(xiJ))return false;
	double xIJ = x[ (j + 1)*w + (i + 1) ]; if(!isfinite(xIJ))return false;
	//if (xij <= 9 || xIj <= 9 || xiJ <= 9 || xIJ <= 9) return false;
	g[0] = 0.5 * (xIj - xij + xIJ - xiJ);
	g[1] = 0.5 * (xiJ - xij + xIJ - xIj);
	//g[0] = xIj - xij;
	//g[1] = xiJ - xij;
	return true;
}

static bool local_orientation_by_tensor(double g[2], struct tdip_state *e,
		int i, int j)
{
	int w = e->w;
	int h = e->h;
	float *x = e->x;
	float *t = e->structure_tensor;
	assert(t);
	if (!insideP(w, h, i, j)) return false;
	if (!insideP(w, h, i+1, j+1)) return false;
	double xij = x[ (j + 0)*w + (i + 0) ]; if(!isfinite(xij))return false;
	double xIj = x[ (j + 0)*w + (i + 1) ]; if(!isfinite(xIj))return false;
	double xiJ = x[ (j + 1)*w + (i + 0) ]; if(!isfinite(xiJ))return false;
	double xIJ = x[ (j + 1)*w + (i + 1) ]; if(!isfinite(xIJ))return false;
	int ij = j*w+i;
	g[0] = t[7*ij + 5];
	g[1] = t[7*ij + 6];
	return true;
}


// generic function to traverse a segment between two pixels
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (qx - px > qy - py || px - qx > qy - py) { // horizontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px + j*slope), j+py, e);
		}
	}
}


struct accumulating_state {
	int w, h;
	float v, *acc;
};

static void pixel_acc(int x, int y, void *ee)
{
	struct accumulating_state *e = ee;
	if (insideP(e->w, e->h, x, y)) {
		//fprintf(stderr, "accu %g in %d %d\n", e->v, x, y);
		int idx = e->w * y + x;
		e->acc[idx] += e->v;
	}
}

static void accumulate_line(float *acc, int w, int h,
		int from[2], int to[2], float v)
{
	struct accumulating_state e = {.w = w, .h = h, .v = v, .acc = acc};
	traverse_segment(from[0], from[1], to[0], to[1], pixel_acc, &e);
}

#include "smapa.h"
SMART_PARAMETER(GRADMIN,9)

void tdip_state_compute(struct tdip_state *e)
{
	// get nice variables
	int w = e->w;
	int h = e->h;
	int tside = e->tside;
	double arad = e->aradius;
	local_orientation_t *local_orientation = e->o;

	// set accumulator to 0
	for (int i = 0; i < tside * tside; i++)
		e->transform[i] = 0;

	// traverse input image and apply writing paradigm
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double g[2] = {0, 0};
		if (local_orientation(g, e, i, j))
		{
			if (fabs(g[1]) < 1e-2 ||  fabs(g[0]) < 1e-2)
				continue;
			double gn = hypot(g[0], g[1]);
			if (gn < GRADMIN())
				continue;

			double theta = i * 2 * M_PI / w;
			double l[3] = {-sin(theta), cos(theta), g[0]/g[1]};

			double rec[4] = {-arad, -arad, arad, arad};
			double aa[2], bb[2]; // (bad names)
			if (cut_line_with_rectangle(aa, bb, l, rec, rec+2))
			{
				double alf = (tside - 1) / (2.0 * arad);
				double bet = (tside - 1) / 2.0;
				int iaa[2] = {alf*aa[0] + bet, alf*aa[1] + bet};
				int ibb[2] = {alf*bb[0] + bet, alf*bb[1] + bet};
				accumulate_line(e->transform, tside, tside,
						iaa, ibb, 1+0*gn);
			}
		}
	}
}

//void tdip(float *transform, double arad, int tside, float *dip, int w, int h)
//{
//	// gradient mask image (for debugging purposes)
//	float *tmp = malloc(w * h * sizeof*tmp);
//	for (int i = 0;  i < w*h; i++)
//		tmp[i] = NAN;
//
//	// set accumulator to 0
//	for (int i = 0; i < tside * tside; i++)
//		transform[i] = 0;
//
//	// traverse input image and apply writing paradigm
//	for (int j = 0; j < h; j++)
//	for (int i = 0; i < w; i++)
//	{
//		double g[2] = {0, 0};
//		if (get_positive_gradient(g, dip, w, h, i, j))
//		{
//			if (fabs(g[1]) < 1e-2 ||  fabs(g[0]) < 1e-2)
//				continue;
//			double gn = hypot(g[0], g[1]);
//			if (gn < GRADMIN())
//				continue;
//			tmp[w*j+i] = gn;
//
//			double theta = i * 2 * M_PI / w;
//			double l[3] = {-sin(theta), cos(theta), g[0]/g[1]};
//
//#if 1
//			double rec[4] = {-arad, -arad, arad, arad};
//			double aa[2], bb[2]; // (bad names)
//			if (cut_line_with_rectangle(aa, bb, l, rec, rec+2))
//			{
//				double alf = (tside - 1) / (2.0 * arad);
//				double bet = (tside - 1) / 2.0;
//				int iaa[2] = {alf*aa[0] + bet, alf*aa[1] + bet};
//				int ibb[2] = {alf*bb[0] + bet, alf*bb[1] + bet};
//				//fprintf(stderr, "i,j=%d,%d, g=%g,%g l=%g,%g,%g, aa=%d,%d bb=%d,%d\n", i, j, g[0], g[1], l[0], l[1], l[2], iaa[0], iaa[1], ibb[0], ibb[1]);
//				accumulate_line(transform, tside, tside,
//						iaa, ibb, 1+0*gn);
//			}
//#endif
//
//#if 0
//			double p[2] = {l[0]/l[2], l[1]/l[2]};
//			int ip[2] = {500 + 40 * p[0], 500 + 40 * p[1]};
//			//fprintf(stderr, "%g %g\t%d %d\n", p[0], p[1], ip[0], ip[1]);
//			if (insideP(tside, tside, ip[0], ip[1]))
//			{
//				double gn = hypot(g[0], g[1]);
//				transform[tside*ip[1] + ip[0]] += gn;
//			}
//#endif
//		}
//	}
//
//	// debug stuff
//	//iio_save_image_float_vec("/tmp/gradmask.tiff", tmp, w, h, 1);
//	free(tmp);
//}

#include "random.c"

void tdip_state_compute_rand_acc(struct tdip_state *e)
{
	// get nice variables
	int w = e->w;
	int h = e->h;
	int tside = e->tside;
	int npairs = e->nrandom;
	double arad = e->aradius;
	local_orientation_t *local_orientation = e->o;

	//(assume the accumulator is already initialized)

	int cx = 0;
	// randomized traversal (independent pairs)
	for (int k = 0; k < npairs; k++)
	{
		int p[2], q[2];
		p[0] = randombounds(0, w-1);
		p[1] = randombounds(0, h-1);
		q[0] = randombounds(0, w-1);
		q[1] = randombounds(0, h-1);
		double vp = e->x[p[1]*w+p[0]];
		double vq = e->x[q[1]*w+q[0]];
		//if (fabs(vp - vq) > 5) continue;
		double gp[2], gq[2];
		if (local_orientation(gp, e, p[0], p[1])
			&& local_orientation(gq, e, q[0], q[1]))
		{
			double ngp = hypot(gp[0], gp[1]);
			double ngq = hypot(gq[0], gq[1]);
			//if (fabs(ngp - ngq) > 2) continue;
			if (fabs(gp[1]) < 1e-2 ||  fabs(gp[0]) < 1e-2)
				continue;
			if (fabs(gq[1]) < 1e-2 ||  fabs(gq[0]) < 1e-2)
				continue;
			double gw = ngp * ngq; // multiplicative weight

			double theta_p = p[0] * 2 * M_PI / w;
			double theta_q = q[0] * 2 * M_PI / w;
			double lp[3]={-sin(theta_p), cos(theta_p), gp[0]/gp[1]};
			double lq[3]={-sin(theta_q), cos(theta_q), gq[0]/gq[1]};

			double x[2];
			cut_two_lines(x, lp, lq);
			double alf = (tside - 1) / (2.0 * arad);
			double bet = (tside - 1) / 2.0;
			int ix[2] = {
				lrint(alf * x[0] + bet),
				lrint(alf * x[1] + bet)
			};
			if (insideP(tside, tside, ix[0], ix[1]))
				e->transform[tside*ix[1] + ix[0]] += gw;
		}
		cx += 1;
	}
	fprintf(stderr, "cx = %d\n", cx);
}


//// accumulates some random samplings into the "transform" image
//void tdipr_acc(float *transform, double arad, int tside, float *dip,
//		int w, int h,
//		int npairs)
//{
//	int cx = 0;
//	// randomized traversal (independent pairs)
//	for (int k = 0; k < npairs; k++)
//	{
//		int p[2], q[2];
//		p[0] = randombounds(0, w-1);
//		p[1] = randombounds(0, h-1);
//		q[0] = randombounds(0, w-1);
//		q[1] = p[1];//randombounds(0, h-1);
//		double vp = dip[p[1]*w+p[0]];
//		double vq = dip[q[1]*w+q[0]];
//		if (fabs(vp - vq) > 5) continue;
//		double gp[2], gq[2];
//		if (get_positive_gradient(gp, dip, w, h, p[0], p[1])
//			&& get_positive_gradient(gq, dip, w, h, q[0], q[1]))
//		{
//			double ngp = hypot(gp[0], gp[1]);
//			double ngq = hypot(gq[0], gq[1]);
//			if (fabs(ngp - ngq) > 2) continue;
//			if (fabs(gp[1]) < 1e-2 ||  fabs(gp[0]) < 1e-2)
//				continue;
//			if (fabs(gq[1]) < 1e-2 ||  fabs(gq[0]) < 1e-2)
//				continue;
//			double gw = ngp * ngq; // multiplicative weight
//
//			double theta_p = p[0] * 2 * M_PI / w;
//			double theta_q = q[0] * 2 * M_PI / w;
//			double lp[3]={-sin(theta_p), cos(theta_p), gp[0]/gp[1]};
//			double lq[3]={-sin(theta_q), cos(theta_q), gq[0]/gq[1]};
//
//			double x[2];
//			cut_two_lines(x, lp, lq);
//			double alf = (tside - 1) / (2.0 * arad);
//			double bet = (tside - 1) / 2.0;
//			int ix[2] = {
//				lrint(alf * x[0] + bet),
//				lrint(alf * x[1] + bet)
//			};
//			if (insideP(tside, tside, ix[0], ix[1]))
//				transform[tside*ix[1] + ix[0]] += 1;//gw;
//		}
//		cx += 1;
//	}
//	fprintf(stderr, "cx = %d\n", cx);
//}

void tdip_state_init(struct tdip_state *e, float *x, int w, int h,
		int tside, double aradius, int nrand, int ntensor)
{
	e->w = w;
	e->h = h;
	e->ntensor = ntensor;
	e->aradius = aradius;
	e->tside = tside;
	e->x = x;
	if (e->ntensor > 0) {
		fprintf(stderr, "ntensor = %d\n", e->ntensor);
		e->o = local_orientation_by_tensor;
		e->structure_tensor = malloc(7 * w * h * sizeof(float));
		compute_structure_tensor_field_ultra_fancy(e->structure_tensor,
				e->x, e->w, e->h, e->ntensor, 100);
		iio_save_image_float_vec("/tmp/debustr.tiff",
				e->structure_tensor, e->w, e->h, 7);
	} else {
		e->o = local_orientation_by_gradient;
		e->structure_tensor = 0;
	}
	e->transform = malloc(e->tside * e->tside * sizeof*e->transform);
	for (int i = 0; i < e->tside * e->tside; i++)
		e->transform[i] = 0;
}

#ifndef OMIT_MAIN_TDIP
#define MAIN_TDIP
#endif//OMIT_MAIN_TDIP

#ifdef MAIN_TDIP
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	bool randomized = pick_option(&c, &v, "r", NULL);
	int nrand = atoi(pick_option(&c, &v, "n", "1000"));
	int tensor = atoi(pick_option(&c, &v, "t", "0"));
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
				"%s aradius arbins <dip >tdip\n", *v);
		//               0  1       2
		return 1;
	}
	double aradius = atof(v[1]);
	int arbins = atoi(v[2]);

	// read input gradient
	int w, h, pd;
	float *dip = iio_read_image_float_vec("-", &w, &h, &pd);
	if (pd != 1) return fprintf(stderr, "I expect a 1-dimensional dip!\n");

	// compute transform
	int tside = 1 + 2 * arbins;

	struct tdip_state e[1];
	tdip_state_init(e, dip, w, h, tside, aradius,
			randomized?nrand:0, tensor);

	if (!randomized)
		tdip_state_compute(e);
	else
		tdip_state_compute_rand_acc(e);
	//tdip_state_recompute(tside, aradius, nrand, tensor);

	//float *transform = malloc(tside * tside * sizeof*transform);
	//tdip_fill_fancy(transform, aradius, tside, dip, w, h, nrand
	//for (int i = 0; i < tside*tside; i++)
	//	transform[i] = 0;

	//if (randomized)
	//	tdipr_acc(transform, aradius, tside, dip, w, h, nrand);
	//else
	//	tdip(transform, aradius, tside, dip, w, h);

	// save output image
	iio_save_image_float_vec("-", e->transform, e->tside, e->tside, 1);

	// cleanup and exti
	free(dip);
	//tdip_state_cleanup(e);
	return 0;
}
#endif//MAIN_TDIP
