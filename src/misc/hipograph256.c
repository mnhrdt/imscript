
void hypograph(float *y, float *x, int w, int h, int n)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int k = 0; k < n; k++)
		y[(j*w+i)*n + k] = x[j*w+i] > k;
}

#include <stdlib.h>
#include "iio.h"
int main(void)
{
	int w, h;
	float *x = iio_read_image_float("-", &w, &h);
	int n = 256;
	float *y = malloc(n*w*h*sizeof*y);
	hypograph(y, x, w, h, n);
	iio_write_image_float_vec("-", y, w, h, n);
	return 0;
}
