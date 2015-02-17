// rpc warper with tiff tile cache

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include <stdio.h>


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

static int tiff_cache_interpolate_float(float *result,
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


void rpc_warpabt(float *outa, float *outb, float *h0, int w, int h, int pd,
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
	eval_rpc(c, rpca, axyh[0], axyh[1], *h0);

	// TODO: compute the stem more inteligently here
	// step = nominal resolution
	double csouth[3];
	eval_rpc(csouth, rpca, axyh[0], axyh[1] + 1, *h0);



	fprintf(stderr, "(%g %g %g) => %.8lf %.8lf\n",
			axyh[0], axyh[1], axyh[2], c[0], c[1]);
	fprintf(stderr, "(%g %g %g) => %.8lf %.8lf\n",
			axyh[0], axyh[1]+1, axyh[2], csouth[0], csouth[1]);
	//double ceast[3];
	//eval_rpc(ceast, rpca, axyh[0] + 1, axyh[1], axyh[2]);
	//fprintf(stderr, "(%g %g %g) => %g %g\n",
	//		axyh[0]+1, axyh[1], axyh[2], ceast[0], ceast[1]);

	//double lon_step = ceast[0] - c[0]; // in degrees
	double lat_step = csouth[1] - c[1]; // in degrees

	double latitude = c[1] * (M_PI/180); // in radians
	double lonfactor = cos(latitude);


	// actually apply the factor
	double lon_step = - lat_step / lonfactor;

	double lon_step_m = (lon_step*lonfactor) * (M_PI/180) * EARTH_RADIUS;
	double lat_step_m = lat_step * (M_PI/180) * EARTH_RADIUS;
	fprintf(stderr, "lon_step = %g (%g meters)\n", lon_step, lon_step_m);
	fprintf(stderr, "lat_step = %g (%g meters)\n", lat_step, lat_step_m);

	if (1) { // some consistency tests
		double pc[3];
		eval_rpci(pc, rpca, c[0], c[1], *h0);
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
		float h0ij = h0[j*w+i];
		double paij[3]; eval_rpci(paij, rpca, lon, lat, h0ij);
		double pbij[3]; eval_rpci(pbij, rpcb, lon, lat, h0ij);
		float *oaij = outa + (j*w + i) * pd;
		float *obij = outb + (j*w + i) * pd;
		tiff_cache_interpolate_float(oaij, ta, paij[0], paij[1]);
		tiff_cache_interpolate_float(obij, tb, pbij[0], pbij[1]);
	}
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
	if (c != 10) {
		fprintf(stderr, "usage:\n\t"
		"%s a.{tiff,rpc} b.{tiff,rpc} ax ay h0.tif o{a,b}.tiff\n", *v);
		//0 1       2    3       4    5  6  7      8   9
		return 1;
	}
	char *filename_a    = v[1];
	char *filename_rpca = v[2];
	char *filename_b    = v[3];
	char *filename_rpcb = v[4];
	double axyh[3] ={atof(v[5]), atof(v[6]), 0};
	char *filename_h0   = v[7];
	char *filename_outa = v[8];
	char *filename_outb = v[9];

	// read input images and rpcs
	//int wa, wb, ha, hb, pd, pdb;
	//float *a  = iio_read_image_float_vec(filename_a, &wa, &ha, &pd);
	//float *b  = iio_read_image_float_vec(filename_b, &wb, &hb, &pdb);
	int megabytes = 800;
	struct tiff_tile_cache ta[1], tb[1];
	tiff_tile_cache_init(ta, filename_a, megabytes);
	tiff_tile_cache_init(tb, filename_b, megabytes);
	int pd = ta->i->spp;
	if (pd != tb->i->spp)
		fail("image color depth mismatch\n");
	
	struct rpc rpca[1];
	struct rpc rpcb[1];
	read_rpc_file_xml(rpca, filename_rpca);
	read_rpc_file_xml(rpcb, filename_rpcb);

	int w, h;
	float *h0 = iio_read_image_float(filename_h0, &w, &h);


	// allocate space for output images
	float *outa = xmalloc(w * h * pd * sizeof*outa);
	float *outb = xmalloc(w * h * pd * sizeof*outb);

	// run the algorithm
	rpc_warpabt(outa, outb, h0, w,h,pd, ta,rpca, tb,rpcb, axyh);

	// save the output images
	iio_save_image_float_vec(filename_outa, outa, w, h, pd);
	iio_save_image_float_vec(filename_outb, outb, w, h, pd);

	// cleanup and exit
	free(outa);
	free(outb);
	tiff_tile_cache_free(ta);
	tiff_tile_cache_free(tb);
	return 0;
}
int main(int c,char*v[]){return main_rpc_warpabt(c,v);}
#endif//MAIN_MNEHS
