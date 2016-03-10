#include <math.h>
#include <stdlib.h>
#include <string.h>

// construct the symmetric boundary of an image
// (assumes the color channels are split)
static void construct_symmetric_boundary(float *xx, int w, int h, int pd)
{
	float (*x)[h][w] = (void*)xx;

	// 4 corners
	for (int l = 0; l < pd; l++)
	{
		float v = x[l][0][0]+x[l][h-1][0]+x[l][0][w-1]+x[l][h-1][w-1];
		x[l][0][0] -= v / 4;
		x[l][h-1][0] -= v / 4;
		x[l][0][w-1] -= v / 4;
		x[l][h-1][w-1] -= v/4;
	}

	// vertical sides
	for (int l = 0; l < pd; l++)
	for (int j = 1; j < h-1; j++)
	{
		float v = x[l][j][0] + x[l][j][w-1];
		x[l][j][0]   -= v / 2;
		x[l][j][w-1] -= v / 2;
	}

	// horizontal sides
	for (int l = 0; l < pd; l++)
	for (int i = 1; i < w-1; i++)
	{
		float v = x[l][0][i] + x[l][h-1][i];
		x[l][0][i]   -= v / 2;
		x[l][h-1][i] -= v / 2;
	}

	// interior
	for (int l = 0; l < pd; l++)
	for (int j = 1; j < h - 1; j++)
	for (int i = 1; i < w - 1; i++)
		x[l][j][i] = NAN;
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel(float *x, int w, int h, int i, int j)
{
	//if (i < 0) i = 0; //(never happens in this program)
	//if (j < 0) j = 0; //(never happens in this program)
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = getpixel(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

// evaluate a bilinear cell at the given point
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel(x, w, h, ip  , iq  );
	float b = getpixel(x, w, h, ip+1, iq  );
	float c = getpixel(x, w, h, ip  , iq+1);
	float d = getpixel(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = bilinear_interpolation(in, iw, ih,
						(i-0.5)/2, (j-0.5)/2);
}

// inpaint the NAN values of an image
void simplest_inpainting(float *x, int w, int h)
{
	float *init = malloc(w * h * sizeof*init);
	if (w > 1 || h > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *xs  = malloc(ws * hs * sizeof*xs);
		zoom_out_by_factor_two(xs, ws, hs, x, w, h);
		simplest_inpainting(xs, ws, hs);
		zoom_in_by_factor_two(init, w, h, xs, ws, hs);
		free(xs);
	}
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			x[i] = init[i];
	free(init);
}

// inpaint the NAN values of an image (with split channels)
void simplest_inpainting_separate(float *in, int w, int h, int pd)
{
	for (int l = 0; l < pd; l++)
		simplest_inpainting(in + w*h*l, w, h);
}

// fill-in the periodic component of an image (with split channels)
void ppsmooth(float *y, float *x, int w, int h, int pd)
{
	memcpy(y, x, w*h*pd*sizeof*x);
	construct_symmetric_boundary(y, w, h, pd);
	simplest_inpainting_separate(y, w, h, pd);
	for (int i = 0; i < w*h*pd; i++)
		y[i] = x[i] - y[i];
}


#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	if ((c != 1 && c != 2 && c != 3) || (c>1 && !strcmp(v[1], "-h"))) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_i = c > 1 ? v[1] : "-";
	char *filename_o = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_i, &w, &h, &pd);
	float *y = malloc(w*h*pd*sizeof*y);

	ppsmooth(y, x, w, h, pd);

	iio_save_image_float_split(filename_o, y, w, h, pd);

	return 0;
}
