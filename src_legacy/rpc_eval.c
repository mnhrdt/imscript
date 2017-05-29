#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"

int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s rpc {1|-1} x y h\n", *v);
		//                          0 1    2     3 4 5
		return 1;
	}
	struct rpc r[1]; read_rpc_file_xml(r, v[1]);
	//print_rpc(stderr, r, "a");
	int direction = atoi(v[2]);
	double x[3] = {atof(v[3]), atof(v[4]), atof(v[5])}, y[2];
	if (1 == direction) {
		eval_rpc(y, r, x[0], x[1], x[2]);
	} else if (-1 == direction) {
		eval_rpci(y, r, x[0], x[1], x[2]);
	} else return fprintf(stderr, "unrecognized direction %d\n", direction);
	printf("%g %g %g\n", y[0], y[1], x[2]);
	return 0;
}
