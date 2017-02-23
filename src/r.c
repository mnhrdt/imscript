#include <stdio.h>
#include "iio.h"

int main(int c, char *v[])
{
	if (c != 2) return 1;
	int w, h;
	float *x = iio_read_image_float(v[1], &w, &h);
	printf("w h x = %d %d %p\n", w, h, (void*)x);
	iio_write_image_float("/tmp/ooo.tiff", x, w, h);
	return 0;
}
