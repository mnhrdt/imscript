#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"

void fill_cylpoint(double out[2], struct rpc *ra, struct rpc *rb, double in[2])
{
	double a[2], b[2];
	eval_rpc_pair(a, ra, rb, in[0], in[1], ra->offset[2] + ra->scale[2]/2);
	eval_rpc_pair(b, ra, rb, in[0], in[1], ra->offset[2] - ra->scale[2]/2);
	double alpha = a[1] - b[1];
	double beta  = b[0] - a[0];
	double gamma = a[0]*b[1] - a[1]*b[0];
	double Cx = (rb->dmval[0] + rb->dmval[2])/2;
	double Cy = (rb->dmval[1] + rb->dmval[3])/2;
	double P[3] = {Cx, Cy, hypot(Cx, Cy)/sqrt(2)};
	double n = hypot(alpha,beta);
	double p = alpha / n;
	double q = beta  / n;
	double r = (alpha*P[0] + P[1] + gamma)/(-P[2]*n);
	out[0] = atan2(p, q);
	out[1] = r;
}


int main(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s rpca rpcb ssf > cylpoints\n", *v);
		//                          0 1    2    3
		return 1;
	}
	struct rpc ra[1]; read_rpc_file_xml(ra, v[1]);
	struct rpc rb[1]; read_rpc_file_xml(rb, v[2]);
	int f = atoi(v[3]);
	int nx = (ra->dmval[2] - ra->dmval[0])/f;
	int ny = (ra->dmval[3] - ra->dmval[1])/f;
	for (int j = 0; j <= ny; j++)
	for (int i = 0; i <= nx; i++)
	{
		double p[2] = {ra->dmval[0] + f*i, ra->dmval[1] + f*j}, q[2];
		fill_cylpoint(q, ra, rb, p);
		printf("%g %g %g %g\n", p[0], p[1], q[0], q[1]);
	}
	return 0;
}
