#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"

// 0 1
// 2 3
static float spectral_norm(float A[4])
{
	float a=A[0], b=A[1], c=A[2], d=A[3];
	// p q    a c  a b
	// q r    b d  c d
	float p = a*a + c*c;
	float q = a*b + c*d;
	float r = b*b + d*d;
	float T = p + r;
	float D = p*r - q*q;
	// t=x+y
	// d=x*y
	// y=t-x
	// d=x*t-x*x
	// x*x-t*x+d=0
	// x=(t+sqrt(t*t-4*d))/2
	float s = (T + sqrt(T*T-4*D))/2;
	return s;
}

// 0 1
// 2 3
static void invert_matrix(float invA[4], float A[4])
{
	float a=A[0], b=A[1], c=A[2], d=A[3];
	float det = a*d - b*c;
	invA[0] = d/det;
	invA[1] = -b/det;
	invA[2] = -c/det;
	invA[3] = a/det;
}

static void flowjac(float *y, float *flow, int w, int h)
{
	float (*out)[w][2] = (void*)y;
	getsample_operator p = getsample_1;
#define U(i,j) p(flow,w,h,2,(i),(j),0)
#define V(i,j) p(flow,w,h,2,(i),(j),1)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		//float ux = U(i,j) - U(i-1,j);
		//float vx = V(i,j) - V(i-1,j);
		//float uy = U(i,j) - U(i,j-1);
		//float vy = V(i,j) - V(i,j-1);
		float ux = (-U(i-1,j-1)-2*U(i-1,j)-U(i-1,j+1)
			   +U(i+1,j-1)+2*U(i+1,j)+U(i+1,j+1))/8;
		float vx = (-V(i-1,j-1)-2*V(i-1,j)-V(i-1,j+1)
			   +V(i+1,j-1)+2*V(i+1,j)+V(i+1,j+1))/8;
		float uy = (-U(i-1,j-1)-2*U(i,j-1)-U(i+1,j-1)
			   +U(i-1,j+1)+2*U(i,j+1)+U(i+1,j+1))/8;
		float vy = (-V(i-1,j-1)-2*V(i,j-1)-V(i+1,j-1)
			   +V(i-1,j+1)+2*V(i,j+1)+V(i+1,j+1))/8;
		float du[4] = {0+ux, uy, vx, 0+vy}, idu[4];
		invert_matrix(idu, du);
		out[j][i][0] = spectral_norm(du);
		out[j][i][1] = spectral_norm(idu);
	}
#undef U
#undef V
}

int main(int c, char *v[])
{
	if (c != 1 && c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return EXIT_FAILURE;
	}
	char *infile = c > 1 ? v[1] : "-";
	char *outfile = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) error("2D vector field expected");
	float *y = xmalloc(2*w*h*sizeof*y);
	flowjac(y, x, w, h);
	iio_write_image_float_vec(outfile, y, w, h, 2);
	free(x);
	return EXIT_SUCCESS;
}
