#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define xmalloc malloc

int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s idx [in [out]]\n", *v);
		//                         0  1    2   3
		return EXIT_FAILURE;
	}
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";
	int w, h, pd;
	void *data = iio_read_image_float_vec(in, &w, &h, &pd);
	float (*x)[pd] = data;
	int component = atoi(v[1]);
	if (component < 0 || component >= pd) {
		fprintf(stderr, "ERROR: " "can not get %dth component of "
				"%dD pixels\n", component, pd);
		return EXIT_FAILURE;
	}
	float *y = xmalloc(w*h*sizeof*y);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i][component];
	iio_save_image_float(out, y, w, h);
	return EXIT_SUCCESS;
}
