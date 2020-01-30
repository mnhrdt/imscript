// cut the first quarter of an image
#include "iio.h"
int main(void)
{
	int w, h, pd;
	float *x = iio_read_image_float_vec("-", &w, &h, &pd);
	for (int j = 0; j < h/2; j++)
	for (int i = 0; i < w/2; i++)
	for (int l = 0; l < pd; l++)
		x[pd*(j*(w/2)+i)+l] = x[pd*(j*w+i)+l];
	iio_write_image_float_vec("-", x, w/2, h/2, pd);
	return 0;
}
