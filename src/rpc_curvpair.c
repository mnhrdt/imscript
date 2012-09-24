#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"


#include "xmalloc.c"

#include "iio.h"

#include "smapa.h"
SMART_PARAMETER(HINC,1)

#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif

static double rpc_hcurv(struct rpc *ra, struct rpc *rb,
		double px, double py, double h)
{
	double q[3][2];
	eval_rpc_pair(q[0], ra, rb, px, py, h - HINC());
	eval_rpc_pair(q[1], ra, rb, px, py, h);
	eval_rpc_pair(q[2], ra, rb, px, py, h + HINC());
	double u[2]; FORI(2) u[i] = q[0][i] - q[1][i];
	double v[2]; FORI(2) v[i] = q[2][i] - q[1][i];
	double w[2]; FORI(2) w[i] = q[2][i] - q[0][i];
	double uv = u[0]*v[0] + u[1]*v[1];
	double nu = hypot(u[0], u[1]);
	double nv = hypot(v[0], v[1]);
	double nw = hypot(w[0], w[1]);
	double x = uv/(nu*nv);
	double sintheta = sqrt(1-x*x);
	double r = nw / (2*sintheta);
	//fprintf(stderr, "c(%g %g) = %g\n", px, py, r);
	return r;
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
	float (*e)[nx][1] = xmalloc(1*nx*ny*sizeof(float));
	for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	{
		double fij[2] = {ra->dmval[0] + f*i, ra->dmval[1] + f*j};
		e[j][i][0] = rpc_hcurv(ra, rb, fij[0], fij[1], h);

	}
	iio_save_image_float_vec("-", **e, nx, ny, 1);
	return 0;
}
