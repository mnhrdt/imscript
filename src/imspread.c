#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int symmetrize_index_inside(int i, int m)
{
	assert( i >= 0 && i < m);
	int r = 0;
	if (i >= m/2) r = i-m;
	if (i < m/2) r = i;
	return r;
}

static float ijradius(int w, int h, int i, int j)
{
	int I = symmetrize_index_inside(i, w);
	int J = symmetrize_index_inside(j, h);
	return hypot(I, J);
}

float spread(float *x, int w, int h)
{
	float allowed_radius = fmin(w,h);
	float allowed_value = 0;//1e-7;
	long double sum = 0;
	long double s = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float r = ijradius(w, h, i, j);
		float v = x[j*w+i];
		if (v > allowed_value && r < allowed_radius
			//	&& fabs(i) > 3 && fabs(j) > 3
				)
		{
			sum += v;
			s += v * r * r;
		}
	}
	return (sqrt(s / sum))/1.4;
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 1 && c != 2) {
		fprintf(stderr, "usage:\n\t%s [img]\n", *v);
		//                          0  1
		return 1;
	}
	char *in_img = c > 1 ? v[1] : "-";
	int w, h;
	float *x = iio_read_image_float(in_img, &w, &h);
	float s = spread(x, w, h);
	printf("%.25lf\n", s);
	free(x);
	return 0;
}
