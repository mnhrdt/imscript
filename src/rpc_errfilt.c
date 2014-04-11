#include <stdio.h>

#include "parsenumbers.c"
#define DONT_USE_TEST_MAIN
#include "rpc.c"


int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
			"%s rpca rpcb < pairs >errvecs\n", *v);
		//        0 1    2
		return 1;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];

	int n;
	double *p = read_ascii_doubles(stdin, &n);
	n /= 4;

	struct rpc ra[1]; read_rpc_file_xml(ra, filename_rpca);
	struct rpc rb[1]; read_rpc_file_xml(rb, filename_rpcb);

	for (int i = 0; i < n; i++)
	{
		double *x = p + 4*i + 0;
		double *y = p + 4*i + 2;
		double e[2];
		double h = rpc_height2(ra, rb, x[0], x[1], y[0], y[1], e);
		printf("%g\t%g\t%g\t%g     \t%lf\t%lf %lf\n",
				x[0], x[1], y[0], y[1], h, e[0], e[1]);
	}

	return 0;
}
