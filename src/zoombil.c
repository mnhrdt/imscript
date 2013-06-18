#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>


#include "fail.c"
#include "xmalloc.c"
#include "bilinear_interpolation.c"


static void bilinear_interpolant_vec(float *y, int yw, int yh,
		float *x, int xw, int xh, int pd)
{
	float (*yy)[yw][pd] = (void*)y;
	float wfactor = xw/(float)yw;
	float hfactor = xh/(float)yh;
	for (int j = 0; j < yh; j++)
	for (int i = 0; i < yw; i++)
	{
		float p = i*wfactor;
		float q = j*hfactor;
		bilinear_interpolation_vec_at(yy[j][i], x, xw, xh, pd, p, q);
	}
}


static bool innerP(int w, int h, int i, int j)
{
	if (i < 0) return false;
	if (j < 0) return false;
	if (i >= w) return false;
	if (j >= h) return false;
	return true;
}


static void downsav2(float *yy, float *xx, int w, int h, int pd, int m, int n)
{
	int W = w/m;
	int H = h/n;
	float (*x)[w][pd] = (void*)xx;
	float (*y)[W][pd] = (void*)yy;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	for (int l = 0; l < pd; l++)
	{
		float num = 0, sum = 0;
		for (int jj = 0; jj < n; jj++)
		for (int ii = 0; ii < m; ii++)
		if (innerP(w, h, i*m+ii, j*n+jj))
		{
			num += 1;
			sum += x[j*n+jj][i*m+ii][l];
		}
		y[j][i][l] = sum/num;
	}
}

void compute_resize_sizes(int *out_w, int *out_h, int w, int h)
{
	int ow = *out_w;
	int oh = *out_h;
	if (ow == -1 && oh == -1) ow = 400;
	if (oh == -1) oh = (h*ow)/w;
	if (ow == -1) ow = (w*oh)/h;
	if (ow < 0 && oh < 0) {
		// if both input numbers are negative, we compute
		// a zoom-to-fit preserving aspect ratio
		ow = -ow;
		oh = -oh;
		if (oh*w > ow*h)
			oh = (h*ow)/w;
		else
			ow = (w*oh)/h;
	} //else fail("incompatible output sizes %d %d", ow, oh);
	assert(ow > 0);
	assert(oh > 0);
	*out_w = ow;
	*out_h = oh;
}

#include "smapa.h"
SMART_PARAMETER_SILENT(ZOOMBIL_DIRECT,0)

void resize_api_vec(float *y, int yw, int yh, float *x, int xw, int xh, int pd)
{
	//fprintf(stderr, "resize %dx%d -> %dx%d\n", xw, xh, yw, yh);
	if ((yw >= xw && yh >= xh) || ZOOMBIL_DIRECT()>0) { // zoom in
		bilinear_interpolant_vec(y, yw, yh, x, xw, xh, pd);
	} else if (yw < xw && yh < xh) { // zoom out
		int m = floor(xw/(float)yw);
		int n = floor(xh/(float)yh);
		int tw = xw/m;
		int th = xh/n;
		assert(yw <= tw); assert(tw <= xw);
		assert(yh <= th); assert(th <= xh);
		float *t = xmalloc(tw*th*pd*sizeof*t);
		downsav2(t, x, xw, xh, pd, m, n);
		//fprintf(stderr, "w: %d -> %d -> %d\n", xw, tw, yw);
		//fprintf(stderr, "h: %d -> %d -> %d\n", xh, th, yh);
		bilinear_interpolant_vec(y, yw, yh, t, tw, th, pd);
		free(t);
	} else fail("anisotropic resize %dx%d => %dx%d not implemented",
			xw, xh, yw, yh);
}

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s width height in out\n", *v);
		return EXIT_FAILURE;
	}
	int ow = atoi(v[1]);
	int oh = atoi(v[2]);

	int w, h, pd;
	void *x = iio_read_image_float_vec(v[3], &w, &h, &pd);

	compute_resize_sizes(&ow, &oh, w, h);
	float *y = xmalloc(ow*oh*pd*sizeof*y);
	resize_api_vec(y, ow, oh, x, w, h, pd);
	//bilinear_interpolant_vec(y, ow, oh, x, w, h, pd);

	iio_save_image_float_vec(v[4], y, ow, oh, pd);

	free(x);
	free(y);

	fprintf(stderr, "%g %g\n", w/(1.0*ow), h/(1.0*oh));
	return EXIT_SUCCESS;
}
