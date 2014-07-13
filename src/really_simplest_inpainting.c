#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
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
				//m += a[k];
				//cx += 1;
				m = a[k];
				cx = 1;
				break;
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

void simplest_inpainting(float *out, float *in, int w, int h, int scale)
{
	float *init = malloc(w * h * sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins  = malloc(ws * hs * sizeof*ins);
		float *outs = malloc(ws * hs * sizeof*outs);
		zoom_out_by_factor_two(ins, ws, hs, in, w, h);
		simplest_inpainting(outs, ins, ws, hs, scale - 1);
		zoom_in_by_factor_two(init, w, h, outs, ws, hs);
		free(ins);
		free(outs);
	}
	//else {
	//	for (int i = 0 ; i < w*h; i++)
	//		init[i] = 0;
	//}
	for (int i = 0 ; i < w*h; i++)
		out[i] = isfinite(in[i]) ? in[i] : init[i];
	free(init);
}

void simplest_inpainting_separable(float *out, float *in, int w, int h, int pd,
		int nscales)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		simplest_inpainting(outl, inl, w, h, nscales);
	}
}

#include "iio.h"
int main(int argc, char *argv[])
{
	if (argc != 4) {
		fprintf(stderr, "usage:\n\t"
		"%s data.png mask.png out.png\n", *argv);
		//0 1        2        3
		return 1;
	}
	char *filename_in = argv[1];
	char *filename_mask = argv[2];
	char *filename_out = argv[3];

	int nscales = 20;

	int w[2], h[2], pd;
	float *in = iio_read_image_float_split(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	float *out = malloc(*w * *h *pd * sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			for (int l = 0; l < pd; l++)
				in[*w**h*l+i] = NAN;

	simplest_inpainting_separable(out, in, *w, *h, pd, nscales);

	iio_save_image_float_split(filename_out, out, *w, *h, pd);

	return 0;
}
