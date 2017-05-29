// print the pairs of points defined by a vector field

#include <math.h>
#include <stdio.h>
#include "iio.h"

int main(int c, char **v)
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [img.tiff] > pairs.txt]\n", *v);
		//                          0  1         2
		return 0;
	}
	char *filename_in  = c > 1 ? v[1] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	if (pd == 2)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p = x[2*(j*w+i)+0];
		float q = x[2*(j*w+i)+1];
		if (isfinite(p) && isfinite(q))
			printf("%d %d %lf %lf\n", i, j, i+p, j+q);
	}

	return 0;
}
