// Optical flow in the altitude domain, using RPC functions
//
// The algorithm lies at the intersection of graphics, vision and remote sensing

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include <stdio.h>

#include "iio.h"

#include "getpixel.c"
#include "bicubic.c"


#define DONT_USE_TEST_MAIN
#include "rpc.c"

#define TIFFU_OMIT_MAIN
#include "tiffu.c"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EARTH_RADIUS 6378000.0


struct projection_state {
	struct rpc *r;
	int w, h;
	double lon_0, lat_0, lon_d, lat_d;
};

static void project(double out[3], struct projection_state *p, double in[3])
{
	double x = p->lon_0 + in[0] * p->lon_d;
	double y = p->lat_0 + in[1] * p->lat_d;
	eval_rpci(out, p->r, x, y, in[2]);
}

static void build_projection_states(
		struct projection_state *pa, struct projection_state *pb,
		struct rpc *ra, struct rpc *rb,
		double axyh[3],
		int w, int h)
{

	// center in geographic coordinates
	double center[3], csouth[3];
	eval_rpc(center, ra, axyh[0], axyh[1]    , axyh[2]);
	eval_rpc(csouth, ra, axyh[0], axyh[1] + 1, axyh[2]);

	// stepsize given by 1 pixel to the south
	// XXX WARNING WRONG TODO FIXME : assumes North-South oriented image (!)
	double lat_step = csouth[1] - center[1];
	double latitude = center[1] * (M_PI/180);
	double lonfactor = cos(latitude);
	double lon_step = -lat_step / lonfactor;

	double lon_step_m = (lon_step*lonfactor) * (M_PI/180) * EARTH_RADIUS;
	double lat_step_m = lat_step * (M_PI/180) * EARTH_RADIUS;
	fprintf(stderr, "projection center (%g %g %g) => (%.8lf %.8lf)\n",
			axyh[0], axyh[1], axyh[2], center[0], center[1]);
	fprintf(stderr, "lon_step = %g (%g meters)\n", lon_step, lon_step_m);
	fprintf(stderr, "lat_step = %g (%g meters)\n", lat_step, lat_step_m);

	// fill-in the fields
	pa->r = ra;
	pb->r = rb;
	pa->w = pb->w = w;
	pa->h = pb->h = h;
	pa->lon_0 = pb->lon_0 = center[0];
	pa->lat_0 = pb->lat_0 = center[1];
	pa->lon_d = pb->lon_d = lon_step;
	pa->lat_d = pb->lat_d = lat_step;
}

static void vertical_direction(double result[2],
		struct projection_state *P, double ijh[3])
{
	double p0[3], ph[3], ijhp[3] = {ijh[0], ijh[1], ijh[2] + 1};
	project(p0, P, ijh);
	project(ph, P, ijhp);
	result[0] = ph[0] - p0[0];
	result[1] = ph[1] - p0[1];
	//fprintf(stderr, "vertical %p (%g %g %g | %g) = (%g %g)\n", P, ijh[0], ijh[1], ijh[2],ijhp[2], result[0], result[1]);
}

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

static double randu(double a, double b)
{
	return a + random_uniform() * (b - a);
}

#include "smapa.h"
SMART_PARAMETER(PM_WINRADIUS,1)
SMART_PARAMETER(PM_TRIALS,5)
SMART_PARAMETER(PM_NITER,3)
SMART_PARAMETER(PM_MIN,-100)
SMART_PARAMETER(PM_MAX,500)

static float eval_cost_pair(
		struct tiff_tile_cache *ta, int ai, int aj,
		struct tiff_tile_cache *tb, int bi, int bj
		)
{
	int pd = ta->i->spp;

	double r = 0;
	int rad = PM_WINRADIUS();
	for (int j = -rad; j <= rad; j++)
	for (int i = -rad; i <= rad; i++)
	{
		void *pa = tiff_tile_cache_getpixel(ta, ai + i, aj + j);
		void *pb = tiff_tile_cache_getpixel(tb, bi + i, bj + j);
		float fa[pd]; convert_pixel_to_float(fa, ta->i, pa);
		float fb[pd]; convert_pixel_to_float(fb, tb->i, pb);
		for (int l = 0; l < pd; l++)
			r = hypot(r, fa[l] - fb[l]);
	}
	return r;
}

static float eval_cost(
		struct projection_state *PA,
		struct projection_state *PB,
		struct tiff_tile_cache *ta,
		struct tiff_tile_cache *tb,
		int i, int j,
		float h
		)
{
	double ijh[3] = {i, j, h}, paijh[3], pbijh[3];
	project(paijh, PA, ijh);
	project(pbijh, PB, ijh);
	int ai = lrint(paijh[0]);
	int aj = lrint(paijh[1]);
	int bi = lrint(pbijh[0]);
	int bj = lrint(pbijh[1]);
	return eval_cost_pair(ta, ai, aj, tb, bi, bj);
}
		

static void random_search(float *cost, float *out_h, float *init_h,
		int w, int h,
		struct projection_state *PA,
		struct projection_state *PB,
		struct tiff_tile_cache *ta,
		struct tiff_tile_cache *tb)
{
	float min_off = PM_MIN();
	float max_off = PM_MAX();
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j * w + i;
		for (int k = 0; k < PM_TRIALS(); k++)
		{
			float new_h = init_h[idx] + randu(min_off, max_off);
			float new_cost = eval_cost(PA, PB, ta, tb, i, j, new_h);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				out_h[idx] = new_h;
			}
		}
	}
}

static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

static void backward_propagation(float *cost, float *height,
		int w, int h,
		struct projection_state *PA,
		struct projection_state *PB,
		struct tiff_tile_cache *ta,
		struct tiff_tile_cache *tb)
{
	int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1}};
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j * w + i;
		for (int n = 0; n < 4; n++)
		{
			int ii = i + neigs[n][0];
			int jj = j + neigs[n][1];
			if (!insideP(w, h, ii, jj)) continue;
			float hh = height[jj*w+ii];
			float new_cost = eval_cost(PA, PB, ta, tb, ii, jj, hh);
			if (new_cost < 1.2 * cost[idx])
			{
				cost[idx] = new_cost;
				height[idx] = hh;
			}
		}
	}
}

static void backward_propagation2(float *cost, float *height,
		int w, int h,
		struct projection_state *PA,
		struct projection_state *PB,
		struct tiff_tile_cache *ta,
		struct tiff_tile_cache *tb)
{
	int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1}};
	for (int i = w-1; i >= 0; i--)
	for (int j = h-1; j >= 0; j--)
	{
		int idx = j * w + i;
		for (int n = 0; n < 4; n++)
		{
			int ii = i + neigs[n][0];
			int jj = j + neigs[n][1];
			if (!insideP(w, h, ii, jj)) continue;
			float hh = height[jj*w+ii];
			float new_cost = eval_cost(PA, PB, ta, tb, ii, jj, hh);
			if (new_cost < 1.2 * cost[idx])
			{
				cost[idx] = new_cost;
				height[idx] = hh;
			}
		}
	}
}

static void backward_propagation_h(float *cost, float *height,
		int w, int h,
		struct projection_state *PA,
		struct projection_state *PB,
		struct tiff_tile_cache *ta,
		struct tiff_tile_cache *tb)
{
	float hdiff[4] = { -4, -1, 1, 4};
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j * w + i;
		for (int n = 0; n < 4; n++)
		{
			float hh = height[idx] + hdiff[n];
			float new_cost = eval_cost(PA, PB, ta, tb, i, j, hh);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				height[idx] = hh;
			}
		}
	}
}

// RPC Patch Match
void pm_rpc(float *out_h, float *init_h, int ow, int oh,
		struct tiff_tile_cache *ta, struct rpc *ra,
		struct tiff_tile_cache *tb, struct rpc *rb,
		double axyh[3])
{
	struct projection_state PA[1], PB[1];
	build_projection_states(PA, PB, ra, rb, axyh, ow, oh);

	float *cost = xmalloc(ow * oh * sizeof*cost);
	for (int i = 0; i < ow * oh; i++)
		cost[i] = INFINITY;

	int niter = PM_NITER();
	for (int iter = 0; iter < niter; iter++)
	{
		fprintf(stderr, "iteration %d/%d\n", iter+1, niter);
		random_search(cost, out_h, init_h, ow, oh, PA, PB, ta, tb);
		backward_propagation(cost, out_h, ow, oh, PA, PB, ta, tb);
	}

	free(cost);
}


#define MAIN_MNEHS

#ifdef MAIN_MNEHS
#include <stdio.h>
#include "iio.h"
#include "xmalloc.c"


int main_rpc_pm(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// input arguments
	if (c != 9) {
		fprintf(stderr, "usage:\n\t"
		"%s a.{tiff,rpc} b.{tiff,rpc} ax ay in.tif out.tif\n", *v);
		//0 1       2    3       4    5  6  7      8
		return 1;
	}
	char *filename_a    = v[1];
	char *filename_rpca = v[2];
	char *filename_b    = v[3];
	char *filename_rpcb = v[4];
	double axyh[3] ={atof(v[5]), atof(v[6]), 0};
	char *filename_h0   = v[7];
	char *filename_out  = v[8];

	// read input images
	int megabytes = 800;
	struct tiff_tile_cache ta[1], tb[1];
	tiff_tile_cache_init(ta, filename_a, megabytes);
	tiff_tile_cache_init(tb, filename_b, megabytes);
	int pd = ta->i->spp;
	if (pd != tb->i->spp) fail("image color depth mismatch\n");

	// read input rpcs
	struct rpc rpca[1];
	struct rpc rpcb[1];
	read_rpc_file_xml(rpca, filename_rpca);
	read_rpc_file_xml(rpcb, filename_rpcb);

	// read initialized raster
	int w, h;
	float *in_h0 = iio_read_image_float(filename_h0, &w, &h);

	// allocate space for output raster
	float *out_h = xmalloc(w * h * sizeof*out_h);

	// run the algorithm
	pm_rpc(out_h, in_h0, w, h, ta, rpca, tb, rpcb, axyh);

	// save the output raster
	iio_save_image_float(filename_out, out_h, w, h);

	// cleanup and exit
	free(in_h0);
	free(out_h);
	tiff_tile_cache_free(ta);
	tiff_tile_cache_free(tb);
	return 0;
}
int main(int c,char*v[]){return main_rpc_pm(c,v);}
#endif//MAIN_MNEHS
