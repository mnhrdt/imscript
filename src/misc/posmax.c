#include <math.h>
#include <stdio.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2) return 1;
	int w, h;
	float *x = iio_read_image_float(v[1], &w, &h);
	int p = -1, q = -1;
	float max = -INFINITY;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float g = x[j*w+i];
		if (isfinite(g) && g > max)
		{
			max = g;
			p = i;
			q = j;
		}
	}
	printf("%d %d\n", p, q);
	return 0;
}
