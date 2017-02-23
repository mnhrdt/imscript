#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bicubic.c"

static void tilt_and_shear(float *y, int ow, int oh,
		float *x, int w, int h, int pd,
		float tilt, float shear)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float p[2] = {i/tilt - shear*j/tilt, j};
		float *v = y + pd*(ow*j + i);
		bicubic_interpolation_boundary(v, x, w, h, pd, p[0], p[1], 0);
	}
}

#include "fail.c"
#include "xmalloc.c"
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
				"%s tilt shear [in.png [out.png]]\n", *v);
		//               0  1    2      3       4
		return 1;
	}
	float tilt = atof(v[1]);
	float shear = atof(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	//if (tilt <= 0 || shear < 0)
	//	fail("negative tilts or shears are not accepted\n");

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	int ow = ceil(fmax(tilt*w + shear*h, tilt*w));
	int oh = h;
	float *y = xmalloc(ow * oh * pd * sizeof*y);

	tilt_and_shear(y, ow, oh, x, w, h, pd, tilt, shear);

	iio_write_image_float_vec(filename_out, y, ow, oh, pd);

	free(x);
	free(y);
	return 0;
}
