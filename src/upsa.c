#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "iio.h"
#include "fail.c"
#include "marching_squares.c"
#include "marching_interpolation.c"
#include "bicubic.c"

static float getsample(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	if (l >= pd) l = pd - 1;
	return x[(i+j*w)*pd + l];
}

static void setsample(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

static float interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	// return a;
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float fade(float x)
{
	return x * x * (3 - 2 * x);
}

static float fadeinv(float x)
{
	return (1 + cbrt(x) - cbrt(1-x))/2;
}

static float interpolate_bilinear_fadeinv(float a, float b, float c, float d,
					float x, float y)
{
	float fx = fadeinv(x);
	float fy = fadeinv(y);
	return interpolate_bilinear(a, b, c, d, fx, fy);
}

static float interpolate_bilinear_fade(float a, float b, float c, float d,
					float x, float y)
{
	float fx = fade(x);
	float fy = fade(y);
	return interpolate_bilinear(a, b, c, d, fx, fy);
}

static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return interpolate_nearest(a, b, c, d, x, y);
	case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	case -2: return interpolate_bilinear_fade(a, b, c, d, x, y);
	case -3: return interpolate_bilinear_fadeinv(a, b, c, d, x, y);
	default: fail("caca de vaca");
	}
	return -1;
}

static void interpolate_vec(float *out, float *x, int w, int h, int pd,
		float p, float q, int m)
{
	if (m == 3) {
		bicubic_interpolation(out, x, w, h, pd, p, q);
	} else {
		int ip = floor(p);
		int iq = floor(q);
		for (int l = 0; l < pd; l++)
		{
			float a = getsample(x, w, h, pd, ip  , iq  , l);
			float b = getsample(x, w, h, pd, ip  , iq+1, l);
			float c = getsample(x, w, h, pd, ip+1, iq  , l);
			float d = getsample(x, w, h, pd, ip+1, iq+1, l);
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			//fprintf(stderr, "p%g q%g ip%d iq%d a%g b%g c%g d%g l%d v%g\n", p, q, ip, iq, a, b, c, d, l, v);
			out[l] = v;
		}
	}
}

float *zoom(float *x, int w, int h, int pd, int n, int zt,
		int *ow, int *oh)
{
	int W = n*w;// - n;
	int H = n*h;// - n;
	float *y = xmalloc(W*H*pd*sizeof*y), nf = n;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		float tmp[pd];
		interpolate_vec(tmp, x, w, h, pd, i/nf, j/nf, zt);
		for (int l = 0; l < pd; l++)
			setsample(y, W, H, pd, i, j, l, tmp[l]);
	}

	*ow = W;
	*oh = H;
	return y;
}

float *zoom_with_offset(float *x, int w, int h, int pd, int n, int zt,
		int *ow, int *oh, float dx, float dy)
{
	int W = n*w;// - n;
	int H = n*h;// - n;
	float *y = xmalloc(W*H*pd*sizeof*y), nf = n;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		float tmp[pd];
		interpolate_vec(tmp, x, w, h, pd, (i-dx)/nf, (j-dy)/nf, zt);
		for (int l = 0; l < pd; l++)
			setsample(y, W, H, pd, i, j, l, tmp[l]);
	}

	*ow = W;
	*oh = H;
	return y;
}

#include "pickopt.c"
int main_upsa(int c, char *v[])
{
	float off_x = atof(pick_option(&c, &v, "x", "0"));
	float off_y = atof(pick_option(&c, &v, "y", "0"));
	if (c < 3 || c > 5) {
		fprintf(stderr, "usage:\n\t%s zoomf zoomtype [in [out]]\n", *v);
		//                            1     2         3   4
		return 1;
	}

	int zf = atoi(v[1]);
	int zt = atoi(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	int ow, od;
	float *y = zoom_with_offset(x,w,h,pd, zf, zt, &ow,&od, off_x,off_y);

	iio_write_image_float_vec(filename_out, y, ow, od, pd);

	free(x);
	free(y);

	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_upsa(c, v); }
#endif

