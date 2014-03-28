#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "xmalloc.c"

static bool interestP(float  *m, int w, int h, int i, int j, float t)
{
	return m[w*j+i] > t;
}

static void fill_ids(double *ids, float *m, int w, int h, int (*t)[3])
{

}

int ccfilt(double *ids, float *m, int w, int h, float t)
{
	int nm = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		int idx = w*j + i;
		ids[idx] = interestP(m, w, h, i, j, t);
		if (ids[idx])
			nm + 1;
	}
	if (nm) {
		int (*t)[3] = xmalloc(nm * sizeof*t), cx = 0;
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
			if (ids[w*j+i]) {
				t[cx][0] = i;
				t[cx][1] = j;
				t[cx][2] = -42;
				cx += 1;
			}
		assert(cx == nm);
		fill_ids(ids, w, h, t);
		for (int i = 0; i < nm; i++) {
			assert(t[i][2] > 0);
			int idx = w*t[i][1] + t[i][0];
			ids[idx] = t[i][2];
		}
		free(t);
	}
	return nm;
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in_mask [out_ids]]\n", *v);
		//                          0  1        2 
		return 1;
	}
	char *filename_mask = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h;
	float *m = iio_read_image_float(filename_mask, &w, &h);
	double *ids = xmalloc(w * h * sizeof*ids);
	float threshold = nextafterf(0);
	int r = ccfilt(ids, m, w, h, threshold);
	iio_save_image_double(ids, filename_out, w, h);

	fprintf(stderr, "r = %d\n", r);

	free(ids);
	free(m);
	return 0;
}
