#include <stdio.h>

#include "iio.h"

#include "xmalloc.c"

#define DONT_USE_TEST_MAIN
#include "rpc.c"

int main(int c, char *v[])
{
	if (c != 9 && c != 10) {
		fprintf(stderr, "usage:\n\t"
			"%s rpca rpcb a0x a0y b0x b0y x y h\n", *v);
		//        0 1    2    3   4   5   6   7 8 9
		return EXIT_FAILURE;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];
	int offset_a[2] = {atoi(v[3]), atoi(v[4])};
	int offset_b[2] = {atoi(v[5]), atoi(v[6])};
	int x[2] = {atoi(v[7]), atoi(v[8])};
	double hbase = atof(v[9]);

	struct rpc rpca[1]; read_rpc_file_xml(rpca, filename_rpca);
	struct rpc rpcb[1]; read_rpc_file_xml(rpcb, filename_rpcb);

	int w = size_a[0];
	int h = size_a[1];
	float (*f)[w][2] = xmalloc(w * h * 2 * sizeof(float));
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double x[2] = {offset_a[0] + i, offset_a[1] + j}, r[2];
		eval_rpc_pair(r, rpca, rpcb, x[0], x[1], hbase);
		double ox[2] = {r[0] - offset_b[0], r[1] - offset_b[1]};
		f[j][i][0] = ox[0] + offset_a[0] - x[0];
		f[j][i][1] = ox[1] + offset_a[1] - x[1];
		//fprintf(stderr, "i=(%d %d) x=(%g %g) ox=(%g %g) f=(%g %g)\n",
		//		i, j,
		//		x[0], x[1],
		//		ox[0], ox[1],
		//		f[j][i][0], f[j][i][1]);
	}

	iio_save_image_float_vec(filename_flow, f[0][0], w, h, 2);

	return 0;
}
