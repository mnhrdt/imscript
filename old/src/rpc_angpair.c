#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"


#include "xmalloc.c"

#include "iio.h"

#include "smapa.h"
SMART_PARAMETER(HINC,1)

static void rpc_hvector(double v[2],
		struct rpc *ra, struct rpc *rb,
		double px, double py, double h)
{
	double q[2][2];
	eval_rpc_pair(q[0], ra, rb, px, py, h);
	eval_rpc_pair(q[1], ra, rb, px, py, h + HINC());
	v[0] = q[1][0] - q[0][0];
	v[1] = q[1][1] - q[0][1];
}

static double angle(double a[2], double b[2])
{
	double na = hypot(a[0], a[1]);
	double nb = hypot(b[0], b[1]);
	double ab = a[0]*b[0] + a[1]*b[1];
	return acos(ab/(na*nb));
}

// rpca: left rpc xml file
// rpcb: right rpc xml file
// ssf: sub-sampling factor
// h: base height
int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s rpca rpcb  ssf h > err.uv\n", *v);
		//                          0 1    2     3   4
		return 1;
	}
	struct rpc ra[1]; read_rpc_file_xml(ra, v[1]);
	struct rpc rb[1]; read_rpc_file_xml(rb, v[2]);
	int f = atoi(v[3]);
	double h = atof(v[4]);
	int nx = (ra->dmval[2] - ra->dmval[0])/f;
	int ny = (ra->dmval[3] - ra->dmval[1])/f;
	fprintf(stderr, "will build image of size %dx%d\n", nx, ny);
	double pbase[2] = {ra->dmval[0] + f*nx/2, ra->dmval[1] + f*ny/2};
	double vbase[2];
	rpc_hvector(vbase, ra, rb, pbase[0], pbase[1], h);
	float (*e)[nx][3] = xmalloc(3*nx*ny*sizeof(float));
	for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	{
		double fij[2] = {ra->dmval[0] + f*i, ra->dmval[1] + f*j};
		double vij[2];
		rpc_hvector(vij, ra, rb, fij[0], fij[1], h);
		e[j][i][0] = vij[0] - vbase[0];
		e[j][i][1] = vij[1] - vbase[1];
		e[j][i][2] = angle(vbase, vij);

	}
	iio_write_image_float_vec("-", **e, nx, ny, 3);
	return 0;
}
