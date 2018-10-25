#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	int n = c > 1 ? atoi(v[1]) : 0;
	int w, h, pd;
	float *x = iio_read_image_float_vec("-", &w, &h, &pd);
	if (n < 0) n = 0;
	if (n >= h) n = h-1;
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		printf("%g%c", x[(n*w+i)*pd+l], l==pd-1 ? '\n' : '\t');
	return 0;
}
