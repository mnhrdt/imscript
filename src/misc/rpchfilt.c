#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#include "parsenumbers.c"
#define DONT_USE_TEST_MAIN
#include "rpc.c"

#include "smapa.h"
SMART_PARAMETER(RELIEF_EXAGERATION,0.5)

static void mercator(double m[2], double x[2])
{
	double ex = RELIEF_EXAGERATION();
	double lon0 = 0;
	double R = 6378100;
	double deg = M_PI/180;
	m[0] = R * (x[0] - lon0) * deg * ex;
	m[1] = R * log( ( 1 + sin(x[1]*deg) ) / cos(x[1]*deg) ) * ex;
}

int main(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s rpca rpcb a0x a0y b0x b0y < pairs >quartets"
		//        0 1    2    3   4   5   6
			"\n", *v);
		return EXIT_FAILURE;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];
	int offset_a[2] = {atoi(v[3]), atoi(v[4])};
	int offset_b[2] = {atoi(v[5]), atoi(v[6])};

	int n;
	double *p = read_ascii_doubles(stdin, &n);
	n /= 4;

	struct rpc ra[1]; read_rpc_file_xml(ra, filename_rpca);
	struct rpc rb[1]; read_rpc_file_xml(rb, filename_rpcb);

	for (int i = 0; i < n; i++)
	{
		double x[2] = {offset_a[0] + p[4*i+0], offset_a[1] + p[4*i+1]};
		double y[2] = {offset_b[0] + p[4*i+2], offset_b[1] + p[4*i+3]};
		double e, h = rpc_height(ra, rb, x[0], x[1], y[0], y[1], &e);
		double q[2], m[2];
		eval_rpc(q, ra, x[0], x[1], h);
		mercator(m, q);
		printf("%lf %lf %lf %lf %lf %lf\n", q[0],q[1], h, e, m[0],m[1]);
	}

	return EXIT_SUCCESS;
}
