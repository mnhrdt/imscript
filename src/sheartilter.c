#include <stdio.h>

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s prog a.png b.png out.tiff\n", *v);
		//                          0 1    2     3     4
		return 1;
	}
	char *disparity_program = v[1];
	char *filename_a = v[2];
	char *filename_b = v[3];
	char *filename_out = v[4];
	return 0;
}
