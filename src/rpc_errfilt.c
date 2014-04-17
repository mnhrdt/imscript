#include <stdio.h>

#include "parsenumbers.c"
#define DONT_USE_TEST_MAIN
#include "rpc.c"

#include "parsenumbers.c"
#include "pickopt.c"

static void projective_map(double y[2], double H[9], double x[2])
{
	double z = H[6]*x[0] + H[7]*x[1] + H[8];
	double ry0 = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	double ry1 = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
	y[0] = ry0;
	y[1] = ry1;
}

int main(int c, char *v[])
{
	char *Htext = pick_option(&c, &v, "h", "1 0 0  0 1 0  0 0 1");
	if (c != 3) {
		erro: fprintf(stderr, "usage:\n\t"
		"%s rpca rpcb [-h \"h1 ... h9\"] < pairs >errvecs\n", *v);
		//        0 1    2
		return 1;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];
	int nm;
	double *H = alloc_parse_doubles(9, Htext, &nm);
	if (nm != 9) goto erro;

	fprintf(stderr, "H = %g %g %g  %g %g %g  %g %g %g\n",
			H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8]);

	int n;
	double *p = read_ascii_doubles(stdin, &n);
	n /= 4;
	
	for (int i = 0; i < n; i++)
		projective_map(p + 4*i + 2, H, p + 4*i + 2);

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
