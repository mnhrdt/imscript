// input: an image of heights
// output: an image of spherical neighbor counts
// parameters: sphere radius, vertical exaggeration
// optional parameter: a preliminary mask


#include <math.h>

static int innerP(int w, int h, int x, int y)
{
	return x >= 0 && y >= 0 && x < w && y < h;
}

static void sphere_count_neighbors(float *out_count, float *x, int w, int h,
		float sphere_radius, float vertical_ex, float *in_mask)
{
	int offset = ceil(sphere_radius);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float *o = out_count + j*w + i;
		*o = 0;
		if (in_mask && !in_mask[j*w+i]) continue;
		int cx = 0;
		float z0 = x[j*w + i]/vertical_ex;
		for (int dy = -offset; dy <= offset; dy++)
		for (int dx = -offset; dx <= offset; dx++) {
			int ii = i + dx;
			int jj = j + dy;
			if (!innerP(w, h, ii, jj)) continue;

			float r = hypot(dx, dy);
			if (r >= sphere_radius) continue;

			float z = x[jj*w + ii]/vertical_ex;
			if (!isfinite(z)) continue;

			if ( !(in_mask && !in_mask[jj*w+ii])
					&& hypot(r, z - z0) < sphere_radius
			   )
				cx += 1;
		}
		*o = cx;
	}
}



#include "iio.h"

#include "fail.c"
#include "xmalloc.c"

int main(int c, char *v[])
{
	if (c != 5 && c != 6) {
		fprintf(stderr, "usage:\n\t"
			"%s rad vex in_height out_count [in_mask]\n", *v);
		//        0 1   2   3         4          5
		return 0;
	}
	float rad = atof(v[1]);
	float vex = atof(v[2]);
	char *filename_in = v[3];
	char *filename_out = v[4];
	char *filename_mask = v[5];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *m = NULL;
	if (filename_mask) {
		int ww, hh;
		m = iio_read_image_float(filename_mask, &ww, &hh);
		if (w != ww || h != hh)
			fail("input image and mask size mismatch");
	}
	float *o = xmalloc(w*h*sizeof*o);

	sphere_count_neighbors(o, x, w, h, rad, vex, m);

	iio_save_image_float(filename_out, o, w, h);

	free(o);
	free(x);
	if (m)
		free(m);

	return 0;
}
