#include <stdio.h>
#include "iio.h"


int main(int c, char *v[])
{
	if (c != 1 && c != 2) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		return 1;
	}
	char *filename = c > 1 ? v[1] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename, &w, &h, &pd);
	printf("%d %d %d\n", w, h, pd);
	return 0;
}
