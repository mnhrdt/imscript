#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"

static void flowjac(float *y, float *flow, int w, int h)
{
	float (*jac)[w][2] = (void*)y;
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
		float d = ux + vy;
		float r = uy + vx;
		jac[j][i][0] = 1 + d + r;
		jac[j][i][1] = 1 - d + r;
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
	iio_save_image_float_vec(outfile, y, w, h, 2);
	free(x);
	return EXIT_SUCCESS;
}
