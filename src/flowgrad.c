#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"

static void flowgrad(float *y, float *flow, int w, int h)
{
	float (*gradient)[w][4] = (void*)y;
	getsample_operator p = getsample_1;
#define FX(i,j) p(flow,w,h,2,(i),(j),0)
#define FY(i,j) p(flow,w,h,2,(i),(j),1)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		// sobel
		float ux = -FX(i-1,j-1)-2*FX(i-1,j)-FX(i-1,j+1)
			+FX(i+1,j-1)+2*FX(i+1,j)+FX(i+1,j+1);
		float uy = -FX(i-1,j-1)-2*FX(i,j-1)-FX(i+1,j-1)
			+FX(i-1,j+1)+2*FX(i,j+1)+FX(i+1,j+1);
		float vx = -FY(i-1,j-1)-2*FY(i-1,j)-FY(i-1,j+1)
			+FY(i+1,j-1)+2*FY(i+1,j)+FY(i+1,j+1);
		float vy = -FY(i-1,j-1)-2*FY(i,j-1)-FY(i+1,j-1)
			+FY(i-1,j+1)+2*FY(i,j+1)+FY(i+1,j+1);
		gradient[j][i][0] = ux;
		gradient[j][i][1] = uy;
		gradient[j][i][2] = vx;
		gradient[j][i][3] = vy;
	}
#undef FX
#undef FY
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
	float *y = xmalloc(4*w*h*sizeof*y);
	flowgrad(y, x, w, h);
	iio_save_image_float_vec(outfile, y, w, h, 4);
	free(x);
	return EXIT_SUCCESS;
}
