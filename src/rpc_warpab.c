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


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EARTH_RADIUS 6378000.0


void rpc_warpab(float *outa, float *outb, int w, int h, int pd,
		float *a, int wa, int ha, struct rpc *rpca,
		float *b, int wb, int hb, struct rpc *rpcb,
		double axyh[3])
{
	// PA = rpca inverse
	// LA = rpca direct
	// PB = rpcb inverse
	// LB = rpcb direct

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
	lon_step /= lonfactor;

	if (1) { // some consistency tests
		double pc[3];
		eval_rpci(pc, rpca, c[0], c[1], axyh[2]);
		fprintf(stderr, "(%g %g %g) => %g %g\n",
				c[0], c[1], axyh[2], pc[0], pc[1]);
	}


	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double lon = c[0] + i * lon_step;
		double lat = c[1] + j * lat_step;
		double paij[3]; eval_rpci(paij, rpca, lon, lat, axyh[2]);
		double pbij[3]; eval_rpci(pbij, rpcb, lon, lat, axyh[2]);
		int oidx = (j*w + i) * pd;
		float *oaij = outa + oidx;
		float *obij = outb + oidx;
		bicubic_interpolation(oaij, a, wa, ha, pd, paij[0], paij[1]);
		bicubic_interpolation(obij, b, wb, hb, pd, pbij[0], pbij[1]);
	}
}


#define MAIN_MNEHS

#ifdef MAIN_MNEHS
#include <stdio.h>
#include "iio.h"
#include "xmalloc.c"

int main_rpc_warpab(int c, char *v[])
{
	// input arguments
	if (c != 12) {
		fprintf(stderr, "usage:\n\t"
		"%s a.{tiff,rpc} b.{tiff,rpc} ax ay h0 w h o{a,b}.tiff\n", *v);
		//0 1       2    3       4    5  6  7  8 9 10  11
		return 1;
	}
	char *filename_a    = v[1];
	char *filename_rpca = v[2];
	char *filename_b    = v[3];
	char *filename_rpcb = v[4];
	double axyh[3] ={atof(v[5]), atof(v[6]), atof(v[7])};
	int  w =         atoi(v[8]);
	int  h =         atoi(v[9]);
	char *filename_outa = v[10];
	char *filename_outb = v[11];

	// read input images and rpcs
	int wa, wb, ha, hb, pd, pdb;
	float *a  = iio_read_image_float_vec(filename_a, &wa, &ha, &pd);
	float *b  = iio_read_image_float_vec(filename_b, &wb, &hb, &pdb);
	struct rpc rpca[1];
	struct rpc rpcb[1];
	read_rpc_file_xml(rpca, filename_rpca);
	read_rpc_file_xml(rpcb, filename_rpcb);

	if (pd != pdb)
		fail("input pair has different color depth");

	// allocate space for output images
	float *outa = xmalloc(w * h * pd * sizeof*outa);
	float *outb = xmalloc(w * h * pd * sizeof*outb);

	// run the algorithm
	rpc_warpab(outa, outb, w,h,pd, a,wa,ha,rpca, b,wb,hb,rpcb, axyh);

	// save the output images
	iio_save_image_float_vec(filename_outa, outa, w, h, pd);
	iio_save_image_float_vec(filename_outb, outb, w, h, pd);

	// cleanup and exit
	free(outa);
	free(outb);
	free(a);
	free(b);
	return 0;
}
int main(int c,char*v[]){return main_rpc_warpab(c,v);}
#endif//MAIN_MNEHS
