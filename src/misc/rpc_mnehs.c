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

static void tiff_cache_interpolate_float(float *result,
		struct tiff_tile_cache *t,
		float x, float y)
{
	int pd = t->i->spp;

	x -= 1;
	y -= 1;
	int ix = floor(x);
	int iy = floor(y);

	float c[4][4][pd];
	for (int j = 0; j < 4; j++)
	for (int i = 0; i < 4; i++)
	{
		void *p = tiff_tile_cache_getpixel(t, ix + i, iy + j);
		convert_pixel_to_float(c[i][j], t->i, p);
	}

	for (int l = 0; l < pd; l++) {
		float C[4][4];
		for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			C[i][j] = c[i][j][l];
		float r = bicubic_interpolation_cell(C, x - ix, y - iy);
		result[l] = r;
	}
}

static void tiff_cache_gradient_float(float *result,
		struct tiff_tile_cache *t,
		float x, float y)
{
	if (t->i->spp != 1)
		fail("only gray images yet");
	float r0, rx, ry, eps = 0.1;
	tiff_cache_interpolate_float(&r0, t, x      , y      );
	tiff_cache_interpolate_float(&rx, t, x + eps, y      );
	tiff_cache_interpolate_float(&ry, t, x      , y + eps);
	result[0] = (rx - r0) / eps;
	result[1] = (ry - r0) / eps;
}

void rpc_warpabt(float *outa, float *outb, int w, int h, int pd,
		struct tiff_tile_cache *ta, struct rpc *rpca,
		struct tiff_tile_cache *tb, struct rpc *rpcb,
		double axyh[3])
{
	// PA = rpca inverse
	// LA = rpca direct
	// PB = rpcb inverse
	// LB = rpcb direct
	// FALSE! it is actually the opposite (?)

	// center in geographic coordinates
	double c[3];
	eval_rpc(c, rpca, axyh[0], axyh[1], axyh[2]);
	fprintf(stderr, "(%g %g %g) => %g %g\n",
			axyh[0], axyh[1], axyh[2], c[0], c[1]);

	// TODO: compute the stem more inteligently here
	// step = nominal resolution
	double csouth[3];
	eval_rpc(csouth, rpca, axyh[0], axyh[1] + 1, axyh[2]);
	fprintf(stderr, "(%g %g %g) => %g %g\n",
			axyh[0], axyh[1]+1, axyh[2], csouth[0], csouth[1]);
	double ceast[3];
	eval_rpc(ceast, rpca, axyh[0] + 1, axyh[1], axyh[2]);
	fprintf(stderr, "(%g %g %g) => %g %g\n",
			axyh[0]+1, axyh[1], axyh[2], ceast[0], ceast[1]);

	double lon_step = ceast[0] - c[0]; // in degrees
	double lat_step = csouth[1] - c[1]; // in degrees

	double latitude = c[1] * (M_PI/180); // in radians
	double lonfactor = cos(latitude);

	double lon_step_m = lon_step * (M_PI/180) * EARTH_RADIUS / lonfactor;
	double lat_step_m = lat_step * (M_PI/180) * EARTH_RADIUS;

	fprintf(stderr, "lon_step = %g (%g meters)\n", lon_step, lon_step_m);
	fprintf(stderr, "lat_step = %g (%g meters)\n", lat_step, lat_step_m);

	// actually apply the factor
	lon_step = lonfactor;

	if (1) { // some consistency tests
		double pc[3];
		eval_rpci(pc, rpca, c[0], c[1], axyh[2]);
		fprintf(stderr, "(%g %g %g) => %g %g\n",
				c[0], c[1], axyh[2], pc[0], pc[1]);
	}

	assert(pd = ta->i->spp);
	assert(pd = tb->i->spp);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double lon = c[0] + i * lon_step;
		double lat = c[1] + j * lat_step;
		double paij[3]; eval_rpci(paij, rpca, lon, lat, axyh[2]);
		double pbij[3]; eval_rpci(pbij, rpcb, lon, lat, axyh[2]);
		float *oaij = outa + (j*w + i) * pd;
		float *obij = outb + (j*w + i) * pd;
		tiff_cache_interpolate_float(oaij, ta, paij[0], paij[1]);
		tiff_cache_interpolate_float(obij, tb, pbij[0], pbij[1]);
	}
}


SMART_PARAMETER(TAU,0.25)
SMART_PARAMETER(ALPHA,1)
SMART_PARAMETER(NITER,1)
SMART_PARAMETER(NWARPS,1)
SMART_PARAMETER(NSCALES,1)

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


// Modèle Numérique d'Élévation Horn Schunck (cas général avec RPC)
void mnehs_rpc(float *out_h, float *init_h, int ow, int oh,
		struct tiff_tile_cache *ta, struct rpc *ra,
		struct tiff_tile_cache *tb, struct rpc *rb,
		double axyh[3],
		float alpha2, int niter)
{
	// numeric convention: float=values, double=positions
	int pd = ta->i->spp;
	if (pd != 1) fail("only gray images by now");

	// compute projection functions with proper normalization
	struct projection_state PA[1], PB[1];
	build_projection_states(PA, PB, ra, rb, axyh, ow, oh);

	// allocate temporary images
	float *h   = xmalloc(ow * oh * sizeof*h);     // h-increment
	float *Q   = xmalloc(ow * oh * sizeof*Q);     // Q
	float *amb = xmalloc(ow * oh * sizeof*amb);   // A-B (warped)
	float *iga = xmalloc(ow * oh * 2 * sizeof*Q); // grad(A) (warped)
	float *igb = xmalloc(ow * oh * 2 * sizeof*Q); // grad(B) (warped)
	
	// fill images q and amb (a minus b)
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float h0 = getsample_nan(init_h, ow, oh, 1, i, j, 0);
		double ijh[3] = {i, j, h0}, paijh[3], pbijh[3];
		project(paijh, PA, ijh);
		project(pbijh, PB, ijh);
		float va[pd], vb[pd], vga[2*pd], vgb[2*pd];
		tiff_cache_interpolate_float(va, ta, paijh[0], paijh[1]);
		tiff_cache_interpolate_float(vb, tb, pbijh[0], pbijh[1]);
		tiff_cache_gradient_float(vga, ta, paijh[0], paijh[1]);
		tiff_cache_gradient_float(vgb, tb, pbijh[0], pbijh[1]);
		double PAp[2]; vertical_direction(PAp, PA, paijh);
		double PBp[2]; vertical_direction(PBp, PB, pbijh);
		float gapa = vga[0] * PAp[0] + vga[1] * PAp[1];
		float gbpb = vgb[0] * PBp[0] + vgb[1] * PBp[1];
		Q  [j*ow+i] = gapa - gbpb;
		amb[j*ow+i] = va[0] - vb[0];
		iga[(j*ow+i)*2+0] = vga[0];
		iga[(j*ow+i)*2+1] = vga[1];
		igb[(j*ow+i)*2+0] = vgb[0];
		igb[(j*ow+i)*2+1] = vgb[1];
	}

	iio_write_image_float_vec("/tmp/mnehs_amb.tiff", amb, ow, oh, 1);
	iio_write_image_float_vec("/tmp/mnehs_Q.tiff", Q, ow, oh, 1);
	iio_write_image_float_vec("/tmp/mnehs_iga.tiff", iga, ow, oh, 2);
	iio_write_image_float_vec("/tmp/mnehs_igb.tiff", igb, ow, oh, 2);

	// save raster grid
	float *raster_grid = xmalloc(2 * ow * oh * sizeof(float));
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		double x = PA->lon_0 + i * PA->lon_d;
		double y = PA->lat_0 + j * PA->lat_d;
		raster_grid[(j*ow+i)*2+0] = x;
		raster_grid[(j*ow+i)*2+1] = y;

	}
	iio_write_image_float_vec("/tmp/mnehs_raster.tiff",raster_grid,ow,oh,2);
	free(raster_grid);

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
			if (fabs(Q[ij]) > 200) // there is boundary
				h[ij] = -amb[ij]/Q[ij];
		}
	}
	iio_write_image_float_vec("/tmp/mnehs_h.tiff", h, ow, oh, 1);

	// update result
	for (int i = 0; i < ow * oh; i++)
		out_h[i] = init_h[i] + h[i];
	
	// cleanup and exit
	free(h);
	free(Q);
	free(amb);
}


#define MAIN_MNEHS

#ifdef MAIN_MNEHS
#include <stdio.h>
#include "iio.h"
#include "xmalloc.c"


int main_rpc_warpabt(int c, char *v[])
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
	float alpha2 = ALPHA()*ALPHA();
	int niter = NITER();
	int nwarps = NWARPS();
	for (int i = 0; i < nwarps; i++)
	{
		mnehs_rpc(out_h, in_h0, w,h,ta,rpca,tb,rpcb,axyh, alpha2,niter);
		memcpy(in_h0, out_h, w*h*sizeof*in_h0);
	}

	// save the output raster
	iio_write_image_float(filename_out, out_h, w, h);

	// cleanup and exit
	free(in_h0);
	free(out_h);
	tiff_tile_cache_free(ta);
	tiff_tile_cache_free(tb);
	return 0;
}
int main(int c,char*v[]){return main_rpc_warpabt(c,v);}
#endif//MAIN_MNEHS
