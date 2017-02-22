#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 6)
		return fprintf(stderr, "usage:\n\t"
				"%s in out_00 out_01 out_10 out_11\n", *v);
	//                        0 1  2      3      4      5
	char *filename_in = v[1];
	char *filename_out[4];
	for (int i = 0; i < 4; i++)
		filename_out[i] = v[2+i];
	int w, h;
	float *o[4], *x = iio_read_image_float(filename_in, &w, &h);
	for (int i = 0; i < 4; i++)
		o[i] = malloc(w*h*sizeof*x/4);
	for (int j = 0; j < h/2; j++)
	for (int i = 0; i < w/2; i++)
	{
		o[0][j*w/2+i] = x[(2*j+0)*w+(2*i+0)];
		o[1][j*w/2+i] = x[(2*j+0)*w+(2*i+1)];
		o[2][j*w/2+i] = x[(2*j+1)*w+(2*i+0)];
		o[3][j*w/2+i] = x[(2*j+1)*w+(2*i+1)];
	}
	for (int i = 0; i < 4; i++)
		iio_save_image_float(filename_out[i], o[i], w/2, h/2);
	return 0;
}
