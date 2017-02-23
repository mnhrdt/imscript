#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"


#include "xmalloc.c"

#include "iio.h"

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
	float (*e)[nx][2] = xmalloc(2*nx*ny*sizeof(float));
	for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	{
		double fij[2] = {ra->dmval[0] + f*i, ra->dmval[1] + f*j};
		double tij[2], rij[2];
		eval_rpc_pair(tij, ra, rb, fij[0], fij[1], h);
		eval_rpc_pair(rij, rb, ra, tij[0], tij[1], h);
		e[j][i][0] = rij[0] - fij[0];
		e[j][i][1] = rij[1] - fij[1];

	}
	iio_write_image_float_vec("-", **e, nx, ny, 2);
	return 0;
}
