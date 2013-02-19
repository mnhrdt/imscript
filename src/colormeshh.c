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

int readMatrix(double H[3][3], char* name) {
  FILE *f = fopen(name,"r"); 

  int r=0;
//  r += fscanf(f, "[");
  r += fscanf(f, "[ %lf %lf %lf ; ", &H[0][0],  &H[0][1], &H[0][2]);
  r += fscanf(f, "%lf %lf %lf ; ", &H[1][0],  &H[1][1], &H[1][2]);
  r += fscanf(f, "%lf %lf %lf ]", &H[2][0],  &H[2][1], &H[2][2]);
//  r += fscanf(f, "]");

   // fail if the number entries is not 9
  if(r!=9) abort();

  fclose(f);
  return 1;
}

static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double z = H[2][0]*x[0] + H[2][1]*x[1] + H[2][2];
	y[0] = (H[0][0]*x[0] + H[0][1]*x[1] + H[0][2]) / z;
	y[1] = (H[1][0]*x[0] + H[1][1]*x[1] + H[1][2]) / z;
}

static double invert_homography(double invH[3][3], double H[3][3])
{
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = (a[4]*a[8]-a[5]*a[7])/det;
	r[1] = (a[2]*a[7]-a[1]*a[8])/det;
	r[2] = (a[1]*a[5]-a[2]*a[4])/det;
	r[3] = (a[5]*a[6]-a[3]*a[8])/det;
	r[4] = (a[0]*a[8]-a[2]*a[6])/det;
	r[5] = (a[2]*a[3]-a[0]*a[5])/det;
	r[6] = (a[3]*a[7]-a[4]*a[6])/det;
	r[7] = (a[1]*a[6]-a[0]*a[7])/det;
	r[8] = (a[0]*a[4]-a[1]*a[3])/det;
	return det;
}

#include "smapa.h"
SMART_PARAMETER_SILENT(IJMESH,0)

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s colors heights rpc Hfile.txt > ply\n", *v);
		//        0 1      2       3   4
		return 1;
	}
	char *fname_colors = v[1];
	char *fname_heights = v[2];
	char *fname_rpc = v[3];
	double H[3][3], invH[3][3];
	readMatrix(H, v[4]);
	invert_homography(invH, H);

	int w, h, pd, ww, hh;
	uint8_t *colors = iio_read_image_uint8_vec(fname_colors, &w, &h, &pd);
	float *heights = iio_read_image_float(fname_heights, &ww, &hh);
	if (w != ww || h != hh) fail("color and height image size mismatch");
	if (pd != 1 && pd != 3) fail("expecting a gray or color image");

	struct rpc r[1]; read_rpc_file_xml(r, fname_rpc);

	printf("ply\n");
	printf("format ascii 1.0\n");
	//printf("format binary_little_endian 1.0\n");
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

	uint8_t (*color)[w][pd] = (void*)colors;
	float (*height)[w] = (void*)heights;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		uint8_t rgb[3];
		for (int k = 0; k < pd; k++) rgb[k] = color[j][i][k];
		for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
		double xy[2] = {i, j}, pq[2];
		apply_homography(pq, invH, xy);
		double xyz[3] = {pq[0], pq[1], height[j][i]};
		if (!IJMESH())
			getxyz(xyz, r, pq[0], pq[1], height[j][i]);
		printf("%.16lf %.16lf %.16lf %d %d %d\n",
				xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
	}
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {j*w+i, (j+1)*w+i, (j+1)*w+i+1, j*w+i+1};
		printf("4 %d %d %d %d\n", q[0], q[1], q[2], q[3]);
	}

	return 0;
}
