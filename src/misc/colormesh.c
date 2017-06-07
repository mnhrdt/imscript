#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


#include "iio.h"

#include "fail.c"

#define DONT_USE_TEST_MAIN
#include "rpc.c"

static void mercator(double m[2], double x[2])
{
	double R = 6378100;
	double deg = M_PI/180;
	m[0] = R * x[0] * deg;
	m[1] = R * log( ( 1 + sin(x[1]*deg) ) / cos(x[1]*deg) );
}

static void getxyz(double xyz[3], struct rpc *r, int i, int j, double h)
{
	double tmp[2];
	eval_rpc(tmp, r, i, j, h);
	mercator(xyz, tmp);
	xyz[2] = h;
}

#include "smapa.h"
SMART_PARAMETER_SILENT(IJMESH,0)
SMART_PARAMETER_SILENT(NOMESH,0)

int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t"
			"%s colors heights rpc offsetx offsety > ply\n", *v);
		//        0 1      2       3   4       5
		return 1;
	}
	char *fname_colors = v[1];
	char *fname_heights = v[2];
	char *fname_rpc = v[3];
	int offsetx = atoi(v[4]);
	int offsety = atoi(v[5]);

	int w, h, pd, ww, hh;
	uint8_t *colors = iio_read_image_uint8_vec(fname_colors, &w, &h, &pd);
	float *heights = iio_read_image_float(fname_heights, &ww, &hh);
	if (w != ww || h != hh) fail("color and height image size mismatch");
	if (pd != 1 && pd != 3) fail("expecting a gray or color image");

	struct rpc r[1]; read_rpc_file_xml(r, fname_rpc);

	int nvert = 0;
	for (int i = 0; i < w*h; i++)
		if (isfinite(heights[i]))
			nvert += 1;
	int omit_mesh = NOMESH() || nvert != w*h;

	printf("ply\n");
	printf("format ascii 1.0\n");
	//printf("format binary_little_endian 1.0\n");
	printf("comment created by cutrecombine\n");
	printf("element vertex %d\n", nvert);
	printf("property float x\n");
	printf("property float y\n");
	printf("property float z\n");
	printf("property uchar red\n");
	printf("property uchar green\n");
	printf("property uchar blue\n");
	if (!omit_mesh) {
		printf("element face %d\n", (w-1)*(h-1));
		printf("property list uchar int vertex_index\n");
	}
	printf("end_header\n");

	uint8_t (*color)[w][pd] = (void*)colors;
	float (*height)[w] = (void*)heights;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(height[j][i])) continue;
		uint8_t rgb[3];
		for (int k = 0; k < pd; k++) rgb[k] = color[j][i][k];
		for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
		double xyz[3] = {j + offsety, i + offsetx, height[j][i]};
		if (!IJMESH())
			getxyz(xyz, r, i + offsetx, j + offsety, height[j][i]);
		printf("%.16lf %.16lf %.16lf %d %d %d\n",
				xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
	}
	if (!omit_mesh)
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {j*w+i, (j+1)*w+i, (j+1)*w+i+1, j*w+i+1};
		printf("4 %d %d %d %d\n", q[0], q[1], q[2], q[3]);
	}

	return 0;
}
