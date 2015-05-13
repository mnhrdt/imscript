#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// compute the norm of a vector
static double vector_norm(double x[3])
{
	return hypot(hypot(x[0], x[1]), x[2]);
}

// compute the distance between two straight lines
// each line is represented by 6 numbers: a point and a direction vector
// if the straight lines are parallel, it returns NAN, instead of their distance
static double distance_between_two_straight_lines(double a[6], double b[6])
{
	double *u = a + 3;
	double *v = b + 3;
	double w[3]; vector_product(w, u, v);
	double nw = vector_norm(w);
	if (!isnormal(nw))
		return NAN;
	double amb[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
	return fabs(scalar_product(w, amb)/nw);
}

static double matrix_inversion(double ia[9], double a[9])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	ia[0] = ( a[4] * a[8] - a[5] * a[7] ) / det;
	ia[1] = ( a[2] * a[7] - a[1] * a[8] ) / det;
	ia[2] = ( a[1] * a[5] - a[2] * a[4] ) / det;
	ia[3] = ( a[5] * a[6] - a[3] * a[8] ) / det;
	ia[4] = ( a[0] * a[8] - a[2] * a[6] ) / det;
	ia[5] = ( a[2] * a[3] - a[0] * a[5] ) / det;
	ia[6] = ( a[3] * a[7] - a[4] * a[6] ) / det;
	ia[7] = ( a[1] * a[6] - a[0] * a[7] ) / det;
	ia[8] = ( a[0] * a[4] - a[1] * a[3] ) / det;
	return det;
}

static double matrix_times_vector(double y[3], double A[9], double x[3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	y[0] = A[0] * x[0]  +  A[1] * x[1]  +  A[2] * x[2];
	y[1] = A[3] * x[0]  +  A[4] * x[1]  +  A[5] * x[2];
	y[2] = A[6] * x[0]  +  A[7] * x[1]  +  A[8] * x[2];
}

static void matrix_product(double ab[9], double a[9], double b[9])
{
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		ab[3*i+j] = 0;
		for (int k = 0; k < 3; k++)
			ab[3*i+j] += a[3*i+k] * b[3*k+j];
	}
}


static void obtain_camera_center(double center[3], double P[12])
{
	// 1. decompose projection matrix P into rotation and translation parts
	//
	// 0 1 2 3
	// 4 5 6 7
	// 8 9 10 11
	double R[9] = {P[0], P[1], P[2], P[4], P[5], P[6], P[8], P[9], P[10]};
	double d[3] = {-P[3], -P[7], -P[11]};

	// 2. compute camera center
	double iR[9];
	matrix_inversion(iR, R);
	matrix_times_vector(center, iR, d);
}

// find the straight line determined by a pixel on the given camera
static void obtain_pixel_direction(double d[3], double P[12], double ij1[3])
{
	// 1. extract the rotation part of P
	//                                 0 1 2 3
	//                                 4 5 6 7
	//                                 8 9 10 11
	double R[9] = {P[0], P[1], P[2], P[4], P[5], P[6], P[8], P[9], P[10]};
	double iR[9]; matrix_inversion(iR, R);

	// 2. compute line of sight of this pixel
	matrix_times_vector(d, iR, ij1);

	// 3. normalize
	double nd = vector_norm(d);
	for (int k = 0; k < 3; k++)
		d[k] /= nd;
}

#include "iio.h"
#include "parsenumbers.c"
#include "pickopt.c"
#include "smapa.h"

SMART_PARAMETER(FRUST_R,255)
SMART_PARAMETER(FRUST_G,0)
SMART_PARAMETER(FRUST_B,0)

int main(int c, char *v[])
{
	// process input arguments
	if (c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s P.txt img.png d_near d_far >out.ply\n", *v);
		//        0 1     2       3      4       5
		return 1;
	}
	char *string_P = v[1];
	char *filename_img = v[2];
	double d_near = atof(v[3]);
	double d_far = atof(v[4]);

	// read input image
	int iw, ih, pd;
	uint8_t *colors = iio_read_image_uint8_vec(filename_img, &iw, &ih, &pd);
	if (pd != 1 && pd != 3)
		return fprintf(stderr, "expecting a gray or color image");

	// read camera matrix
	double P[12];
	read_n_doubles_from_string(P, string_P, 12);

	// setup combinatorial variables
	int w = 100;
	int h = 100;

	int nvertices = 1 + 2 * w * h; // two images, plus the focal point
	int nfaces = 1 + 2 * (w-1) * (h-1); // two images, plus one triangle
	uint8_t (*color)[iw][pd] = (void*)colors;


	// setup geometrical variables
	double center[3];
	obtain_camera_center(center, P);

	// print PLY header
	printf("ply\n");
	printf("format ascii 1.0\n");
	printf("comment created by frustumize\n");
	printf("comment P = %g %g %g %g %g %g %g %g %g %g %g %g\n",
			P[0], P[1], P[2], P[3], P[4], P[5],
			P[6], P[7], P[8], P[9], P[10], P[11]);
	printf("comment frustum interval = %g %g\n", d_near, d_far);
	printf("comment img = %s\n", filename_img);
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

	// setup magic indexes of relevant points
	int idx_focal = 2 * w * h;
	int idx_right = 2 * w * h - 1;
	int idx_left  = idx_right - 2 * w + 2;

	// put the two images on the faces of the frustum
	int idx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		// coordinates in the original image
		int ii = (i * iw) / w;
		int jj = (j * ih) / h;

		// get rgb at this point
		uint8_t rgb[3];
		for (int k = 0; k < pd; k++)
			rgb[k] = color[jj][ii][k];
		for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
		if (idx+1 == idx_left || idx+1 == idx_right)
		{
			rgb[0] = FRUST_R();
			rgb[1] = FRUST_G();
			rgb[2] = FRUST_B();
		}

		// get pixel direction
		double direction[3];
		double ij1[3] = {ii, jj, 1};
		obtain_pixel_direction(direction, P, ij1);

		// cut two sides of the frustum
		double ij_near[3], ij_far[3];
		for (int k = 0; k < 3; k++)
		{
			ij_near[k] = center[k] + d_near * direction[k];
			ij_far[k]  = center[k] + d_far * direction[k];
		}

		// paint the two pixels
		printf("%.16lf %.16lf %.16lf %d %d %d\n",
				ij_near[0], ij_near[1], ij_near[2],
				rgb[0], rgb[1], rgb[2]);
		printf("%.16lf %.16lf %.16lf %d %d %d\n",
				ij_far[0], ij_far[1], ij_far[2],
				rgb[0], rgb[1], rgb[2]);

		idx += 2;
	}
	assert(idx == 2 * w * h);

	// put the focal point (in red, at position "2*w*h")
	printf("%.16lf %.16lf %.16lf %d %d %d\n",
				center[0], center[1], center[2],
				(int)FRUST_R(), (int)FRUST_G(), (int)FRUST_B());


	// base of the triangle
	printf("3 %d %d %d\n", idx_focal, idx_left, idx_right);
	// connect all the pixels
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {j*w+i, (j+1)*w+i, (j+1)*w+i+1, j*w+i+1};
		for (int k = 0; k < 2; k++)
			printf("4 %d %d %d %d\n",
					2*q[0]+k, 2*q[1]+k, 2*q[2]+k, 2*q[3]+k);
	}

	// cleanup and exit
	free(colors);
	return 0;
}
