#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"


#include "xmalloc.c"

#include "iio.h"

// rpc: rpc file in XML
// ssf: sub sampling factor (wrt to the original panchromatic image)
// h: a height at which to evaluate the numbers
int main(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s rpc ssf h > err.uv\n", *v);
		//                          0 1   2   3
		return 1;
	}
	struct rpc r[1]; read_rpc_file_xml(r, v[1]);
	int f = atoi(v[2]);
	double h = atof(v[3]);
	int nx = (r->dmval[2] - r->dmval[0])/f;
	int ny = (r->dmval[3] - r->dmval[1])/f;
	//fprintf(stderr, "
	fprintf(stderr, "will build image of size %dx%d\n", nx, ny);
	float (*e)[nx][2] = xmalloc(2*nx*ny*sizeof(float));
	for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	{
		double fij[2] = {r->dmval[0] + f*i, r->dmval[1] + f*j};
		double tij[2];
		eval_rpc_pair(tij, r, r, fij[0], fij[1], h);
		e[j][i][0] = tij[0] - fij[0];
		e[j][i][1] = tij[1] - fij[1];

	}
	iio_write_image_float_vec("-", **e, nx, ny, 2);
	return 0;
}
