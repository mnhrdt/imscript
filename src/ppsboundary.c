#include <math.h>
#include "iio.h"
int main(void)
{
	int w, h, pd;
	float *xx = iio_read_image_float_vec("-", &w, &h, &pd);
	float (*x)[w][pd] = (void*)xx;

	// 4 corners
	for (int l = 0; l < pd; l++)
	{
		float v = x[0][0][l]+x[h-1][0][l]+x[0][w-1][l]+x[h-1][w-1][l];
		x[0][0][l] -= v / 4;
		x[h-1][0][l] -= v / 4;
		x[0][w-1][l] -= v / 4;
		x[h-1][w-1][l] -= v/4;
	}

	// vertical sides
	for (int j = 1; j < h-1; j++)
	for (int l = 0; l < pd; l++)
	{
		float v = x[j][0][l] + x[j][w-1][l];
		//x[j][0][l] = x[j][w-1][l] = v / 2;
		x[j][0][l]   -= v / 2;
		x[j][w-1][l] -= v / 2;
	}

	// horizontal sides
	for (int i = 1; i < w-1; i++)
	for (int l = 0; l < pd; l++)
	{
		float v = x[0][i][l] + x[h-1][i][l];
		//x[0][i][l] = x[h-1][i][l] = v / 2;
		x[0][i][l]   -= v / 2;
		x[h-1][i][l] -= v / 2;
	}

	// interior
	for (int j = 1; j < h - 1; j++)
	for (int i = 1; i < w - 1; i++)
	for (int l = 0; l < pd; l++)
		x[j][i][l] = NAN;

	iio_save_image_float_vec("-", xx, w, h, pd);
	return 0;
}
