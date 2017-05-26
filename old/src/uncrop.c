#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool uncrop(int *ucx, int *ucy,
		float *aa, int wa, int ha, float *bb, int wb, int hb)
{
	if (wa > wb || ha > hb)
		return false;

	float (*a)[wa] = (void*)aa;
	float (*b)[wb] = (void*)bb;
	for (int i = 0; i < wb-wa; i++)
	for (int j = 0; j < hb-ha; j++) {
		bool matchij = true;
		for (int p = 0; p < wa; p++)
		for (int q = 0; q < ha; q++) {
			assert(i+p < wb);
			assert(j+q < hb);
			if (a[q][p] != b[q+j][p+i]) {
				matchij = false;
				goto endij;
			}
		}
endij:
		if (matchij) {
			*ucx = i;
			*ucy = j;
			return true;
		}
	}
	return false;
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s needle haystack\n", *v);
		return 2;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];

	int wa, ha, wb, hb, ucx, ucy;
	float *a = iio_read_image_float(filename_a, &wa, &ha);
	float *b = iio_read_image_float(filename_b, &wb, &hb);

	bool r = uncrop(&ucx, &ucy, a, wa, ha, b, wb, hb);

	if (r)
		printf("%d %d\n", ucx, ucy);

	return r?0:1;
}
