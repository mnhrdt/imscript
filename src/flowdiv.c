#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"

static void flowdiv(float *y, float *flow, int w, int h)
{
	float (*divergence)[w] = (void*)y;
	getsample_operator p = getsample_1;
#define FX(i,j) p(flow,w,h,2,(i),(j),0)
#define FY(i,j) p(flow,w,h,2,(i),(j),1)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float xx = -FX(i-1,j-1)-2*FX(i-1,j)-FX(i-1,j+1)
			+FX(i+1,j-1)+2*FX(i+1,j)+FX(i+1,j+1);
		float yy = -FY(i-1,j-1)-2*FY(i,j-1)-FY(i+1,j-1)
			+FY(i-1,j+1)+2*FY(i,j+1)+FY(i+1,j+1);
		divergence[j][i] = (xx + yy)/8;
		//divergence[j][i] = xx + yy;
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
	float *y = xmalloc(w*h*sizeof*y);
	flowdiv(y, x, w, h);
	iio_write_image_float_vec(outfile, y, w, h, 1);
	free(x);
	return EXIT_SUCCESS;
}
