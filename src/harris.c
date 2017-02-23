#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"


#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

static void harris(float *yy, float *xx, int w, int h, int pd, float kappa)
{
	float (*x)[w][pd] = (void*)xx;
	float (*y)[w] = (void*)yy;
	FORJ(h) FORI(w) {
		float cornerness = 0;
		FORL(pd) {
			float p =  getsample_0(xx, w, h, pd, i  , j  , l);
			float px = getsample_0(xx, w, h, pd, i+1, j  , l);
			float pmx= getsample_0(xx, w, h, pd, i-1, j  , l);
			float py = getsample_0(xx, w, h, pd, i  , j+1, l);
			float pmy= getsample_0(xx, w, h, pd, i  , j-1, l);
			float lap = 4*p - px - pmx - py - pmy;
			float ctr = hypot(px - pmx, py - pmy);
			//float ix = px - p;
			//float iy = py - p;
			//float H[2][2] = {{ix*ix, ix*iy}, {ix*iy, iy*iy}};
			//float det = H[0][0] * H[1][1] - H[1][0] * H[0][1];
			//// det=0 always !, must consider a window for this to
			//// make sense
			//float tr = H[0][0] + H[1][1];
			//cornerness += det - kappa * tr * tr;
			cornerness += lap - kappa * ctr;
		}
		y[j][i] = cornerness;
	}
}

int main(int c, char *v[])
{
	if (c != 4 && c != 3 && c != 2) {
		fprintf(stderr, "usage:\n\t%s kappa [in [out]]\n", *v);
		//                          0 1  2   3
		return EXIT_FAILURE;
	}
	float kappa = atof(v[1]);
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	float *y = xmalloc(w*h*sizeof*y);
	harris(y, x, w, h, pd, kappa);
	iio_write_image_float(out, y, w, h);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
