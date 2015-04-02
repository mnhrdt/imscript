#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "iio.h"

int main(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s colors heights > ply\n", *v);
		//                          0 1      2
		return 1;
	}
	char *fname_colors = v[1];
	char *fname_heights = v[2];

	// read input images
	int w, h, pd, ww, hh;
	uint8_t *colors = iio_read_image_uint8_vec(fname_colors, &w, &h, &pd);
	float *heights = iio_read_image_float(fname_heights, &ww, &hh);
	if (w != ww || h != hh)
		return fprintf(stderr, "color and height image size mismatch");
	if (pd != 1 && pd != 3)
		return fprintf(stderr, "expecting a gray or color image");

	// assign comfortable pointers
	uint8_t (*color)[w][pd] = (void*)colors;
	float (*height)[w] = (void*)heights;
	int (*vid)[w] = malloc(w*h*sizeof(int));

	// count number of valid vertices
	int nvertices = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isfinite(height[j][i]))
			vid[j][i] = nvertices++;
		else
			vid[j][i] = -1;

	// count number of valid faces
	int nfaces = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 3)
			nfaces += 1;
	}

	// print header
	printf("ply\n");
	printf("format ascii 1.0\n");
	//printf("format binary_little_endian 1.0\n");
	printf("comment created by cutrecombine\n");
	printf("element vertex %d\n", nvertices);
	printf("property float x\n");
	printf("property float y\n");
	printf("property float z\n");
	printf("property uchar red\n");
	printf("property uchar green\n");
	printf("property uchar blue\n");
	printf("element face %d\n", nfaces);
	printf("property list uchar int vertex_index\n");
	printf("end_header\n");
	int cx;

	// output vertices
	cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(height[j][i])) continue;
		uint8_t rgb[3];
		for (int k = 0; k < pd; k++) rgb[k] = color[j][i][k];
		for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
		double xyz[3] = {i, j, height[j][i]};
		printf("%.16lf %.16lf %.16lf %d %d %d\n",
				xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
		cx += 1;
	}
	assert(cx == nvertices);

	cx = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 3)
		{
			printf("4 %d %d %d %d\n", q[0], q[1], q[2], q[3]);
			cx += 1;
		}
	}
	assert(cx == nfaces);

	return 0;
}
