#include <stdio.h>
#include <stdlib.h>

#include "iio.h"

// deinterlace and produce two images
int main(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s in out1 out2\n", *v);
		//                          0 1  2    3
		return EXIT_FAILURE;
	}

	int w, h, pd;
	void *xx = iio_read_image_float_vec(v[1], &w, &h, &pd);
	float (*x)[w][pd] = xx;
	float (*y)[w][pd] = calloc(w*h*pd,sizeof(float));
	float (*z)[w][pd] = calloc(w*h*pd,sizeof(float));

	for (int j = 1; j < h; j += 2) {
		for (int i = 0; i < w; i++)
		for (int l = 0; l < pd; l++)
			z[j][i][l] = x[j][i][l];
		if (j+2 < h)
			for (int i = 0; i < w; i++)
			for (int l = 0; l < pd; l++)
				z[j+1][i][l] = (x[j][i][l] + x[j+2][i][l])/2;
	}
	for (int j = 0; j < h; j += 2) {
		for (int i = 0; i < w; i++)
		for (int l = 0; l < pd; l++)
			y[j][i][l] = x[j][i][l];
		if (j+2 < h)
			for (int i = 0; i < w; i++)
			for (int l = 0; l < pd; l++)
				y[j+1][i][l] = (x[j][i][l] + x[j+2][i][l])/2;
	}

	// the other way to interpet it:
//	for (int j = 0; j < h-1; j += 2) {
//		for (int i = 0; i < w; i++)
//		for (int l = 0; l < pd; l++)
//			z[j][i][l] = x[j+1][i][l];
//		if (j+3 < h)
//			for (int i = 0; i < w; i++)
//			for (int l = 0; l < pd; l++)
//				z[j+1][i][l] = (x[j+1][i][l] + x[j+3][i][l])/2;
//	}
//	for (int j = 0; j < h; j += 2) {
//		for (int i = 0; i < w; i++)
//		for (int l = 0; l < pd; l++)
//			y[j][i][l] = x[j][i][l];
//		if (j+2 < h)
//			for (int i = 0; i < w; i++)
//			for (int l = 0; l < pd; l++)
//				y[j+1][i][l] = (x[j][i][l] + x[j+2][i][l])/2;
//	}

	iio_save_image_float_vec(v[2], y[0][0], w, h, pd);
	iio_save_image_float_vec(v[3], z[0][0], w, h, pd);

	return EXIT_SUCCESS;
}
