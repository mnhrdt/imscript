#include <stdio.h>
#include "iio.h"

int main(void)
{
	int w, h;
	float *x = iio_read_image_float("-", &w, &h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (x[w*j+i] > 0)
			printf("%d %d\n", i, j);
	return 0;
}
