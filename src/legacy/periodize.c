#include "iio.h"

#include "xmalloc.c"

int main(int c, char *v[])
{
	int w, h, pd;
	float *x_raw = iio_read_image_float_vec("-", &w, &h, &pd);
	float *y_raw = xmalloc(4*w*h*pd*sizeof(float));
	float (*x)[w][pd] = (void*)x_raw;
	float (*y)[2*w][pd] = (void*)y_raw;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
	{
		float g = x[j][i][l];
		y[j][i][l] = g;
		y[2*h-j-1][i][l] = g;
		y[j][2*w-i-1][l] = g;
		y[2*h-j-1][2*w-i-1][l] = g;
	}
	iio_write_image_float_vec("-", y_raw, 2*w, 2*h, pd);
	return 0;
}
