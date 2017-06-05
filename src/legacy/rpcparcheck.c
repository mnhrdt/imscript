#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#include "xmalloc.c"
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

static void *list_of_first_points(int *n, double cropa[4], double hinte[2])
{
	int N = 69;
	*n = N*N+4*N;
	//double lp[8][3] = {
	//	{cropa[0], cropa[1], hinte[0]},
	//	{cropa[0], cropa[1], hinte[1]},
	//	{cropa[2], cropa[1], hinte[0]},
	//	{cropa[2], cropa[1], hinte[1]},
	//	{cropa[0], cropa[3], hinte[0]},
	//	{cropa[0], cropa[3], hinte[1]},
	//	{cropa[2], cropa[3], hinte[0]},
	//	{cropa[2], cropa[3], hinte[1]}
	//};
	double (*lp)[3] = xmalloc(*n*3*sizeof(double));
	int cx = 0;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	{
		double x = cropa[0] + i*(cropa[2]-cropa[0])/(N-1);
		double y = cropa[1] + j*(cropa[3]-cropa[1])/(N-1);
		lp[cx][0] = x;
		lp[cx][1] = y;
		lp[cx][2] = hinte[0];
		cx += 1;
	}
	for (int i = 0; i < N; i++) {
		lp[cx][0] = cropa[0]; lp[cx][1] = cropa[1];
		lp[cx][2] = hinte[0] + i*(hinte[1]-hinte[0])/(N-1);
		cx += 1;
	}
	for (int i = 0; i < N; i++) {
		lp[cx][0] = cropa[2]; lp[cx][1] = cropa[1];
		lp[cx][2] = hinte[0] + i*(hinte[1]-hinte[0])/(N-1);
		cx += 1;
	}
	for (int i = 0; i < N; i++) {
		lp[cx][0] = cropa[0]; lp[cx][1] = cropa[3];
		lp[cx][2] = hinte[0] + i*(hinte[1]-hinte[0])/(N-1);
		cx += 1;
	}
	for (int i = 0; i < N; i++) {
		lp[cx][0] = cropa[2]; lp[cx][1] = cropa[3];
		lp[cx][2] = hinte[0] + i*(hinte[1]-hinte[0])/(N-1);
		cx += 1;
	}
	assert(cx == *n);
	return lp;
}

static void construct_straight_pairs(double (*pairs)[4],
		struct rpc *ra, struct rpc *rb,
		double (*points)[3], int n)
{
	for (int i = 0; i < n; i++)
	{
		double *point = points[i];
		double *pair = pairs[i];
		double tmp[2];
		pair[0] = point[0];
		pair[1] = point[1];
		eval_rpc(tmp, ra, point[0], point[1], point[2]);
		eval_rpci(pair+2, rb, tmp[0], tmp[1], point[2]);
		if (1) {
			eval_rpc(tmp, rb, pair[2], pair[3], point[2]);
			double tmp2[2];
			eval_rpci(tmp2, ra, tmp[0], tmp[1], point[2]);
			double ee = hypot(tmp2[0]-pair[0], tmp2[1]-pair[1]);
			fprintf(stderr, "epa[%d] = %lf\n", i, ee);
		}
	}
}

static void apply_eucli(double y[2], double x[2], double eucli[5])
{
	double c[2] = {eucli[0], eucli[1]};
	double alpha = eucli[2]*M_PI/180;
	double t[2] = {eucli[3], eucli[4]};
	double xmc[2] = {x[0] - c[0], x[1] - c[1]};
	double rxmc[2] = {
		xmc[0] * cos(alpha) - xmc[1]*sin(alpha),
		xmc[0] * sin(alpha) + xmc[1]*cos(alpha)
	};
	double rx[2] = {rxmc[0] + c[0], rxmc[1] + c[1]};
	y[0] = rx[0] + t[0];
	y[1] = rx[1] + t[1];
}

static void deform_pairs(double (*pairs)[4], int n, double eucli[5])
{
	for (int i = 0; i < n; i++)
	{
		double *x = 2 + pairs[i];
		//fprintf(stderr, "x[i] before = %lf %lf\n", x[0], x[1]);
		apply_eucli(x, x, eucli);
		//fprintf(stderr, "x[i] after = %lf %lf\n", x[0], x[1]);
	}
}

static void build_3dpoints_from_pairs(double (*cloud)[3],
		struct rpc *ra, struct rpc *rb,
		double (*pairs)[4], int n)
{
	for (int i = 0; i < n; i++)
	{
		double *p = pairs[i];
		double e, h = rpc_height(ra, rb, p[0], p[1], p[2], p[3], &e);
		fprintf(stderr, "e[%d] = %lf\n", i, e);
		double x[2];
		eval_rpc(x, ra, p[0], p[1], h);
		mercator(x, x);
		cloud[i][0] = x[0];
		cloud[i][1] = x[1];
		cloud[i][2] = h;
	}
}

int main(int c, char *v[])
{
	if (c != 14) {
		fprintf(stderr, "usage:\n\t"
	"%s rpca rpcb a0x a0y afx afy h0 hf p_rcx p_rcy p_a p_tx p_ty >3dpoints"
	//0 1    2    3   4   5   6   7  8  9     10    11  12   13
			"\n", *v);
		return EXIT_FAILURE;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];
	double cropa[4] = {atof(v[3]), atof(v[4]), atof(v[5]), atof(v[6])};
	double hinte[2] = {atof(v[7]), atof(v[8])};
	double eucli[5] = {atof(v[9]), atof(v[10]), atof(v[11]),
						atof(v[12]), atof(v[13])};

	struct rpc ra[1]; read_rpc_file_xml(ra, filename_rpca);
	struct rpc rb[1]; read_rpc_file_xml(rb, filename_rpcb);

	int n;
	double (*lp)[3] = list_of_first_points(&n, cropa, hinte);
	double (*pairs)[4] = xmalloc(n*4*sizeof(double));
	construct_straight_pairs(pairs, ra, rb, lp, n);
	deform_pairs(pairs, n, eucli);
	double (*cloud)[3] = xmalloc(n*3*sizeof(double));
	build_3dpoints_from_pairs(cloud, ra, rb, pairs, n);

	for (int i = 0; i < n; i++)
		printf("%lf %lf %lf\n", cloud[i][0], cloud[i][1], cloud[i][2]);

	free(lp);
	free(pairs);
	free(cloud);

	return EXIT_SUCCESS;
}
