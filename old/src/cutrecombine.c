#include <stdio.h>

#include "iio.h"

#include "fail.c"

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s colors heights > ply\n", *v);
		//                          0 1      2
		return 1;
	}
	char *filename_colors = v[1];
	char *filename_heights = v[2];

	int w, h, pd, ww, hh;
	float *colors = iio_read_image_float_vec(filename_colors, &w, &h, &pd);
	float *heights = iio_read_image_float(filename_heights, &ww, &hh);
	if (w != ww || h != hh) fail("color and height image size mismatch");
	if (pd != 3) fail("expecting a color image");

	printf("ply\n");
	printf("format ascii 1.0\n");
	printf("comment created by cutrecombine\n");
	printf("element vertex %d\n", w*h);
	printf("property float x\n");
	printf("property float y\n");
	printf("property float z\n");
	printf("property uchar red\n");
	printf("property uchar green\n");
	printf("property uchar blue\n");
	printf("element face %d\n", (w-1)*(h-1));
	printf("property list uchar int vertex_index\n");
	printf("end_header\n");

	float (*color)[w][pd] = (void*)colors;
	float (*height)[w] = (void*)heights;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		printf("%d %d %g %g %g %g\n", i, -j, height[j][i],
				color[j][i][0], color[j][i][1], color[j][i][2]);
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {j*w+i, (j+1)*w+i, (j+1)*w+i+1, j*w+i+1};
		printf("4 %d %d %d %d\n", q[0], q[1], q[2], q[3]);
	}

	return 0;
}
