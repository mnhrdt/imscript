#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



// utility function that always returns a valid pointer to memory
static void *xmalloc(size_t n)
{
	void *new = malloc(n);
	if (!new)
	{
		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
		exit(1);
	}
	return new;
}


// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
inline static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w];
}

// extrapolate by nearest value
inline static float getpixel_1(float *x, int w, int h, int i, int j)
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
	getpixel_operator p = getpixel_1;

	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i  , 2*j  );
		a[1] = p(in, iw, ih, 2*i+1, 2*j  );
		a[2] = p(in, iw, ih, 2*i  , 2*j+1);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

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

static float bilinear_interpolation(float *img, int w, int h, float x, float y)
{
	getpixel_operator p = getpixel_1;

	int ix = x;
	int iy = y;
	float a = p(img, w, h, ix  , iy  );
	float b = p(img, w, h, ix+1, iy  );
	float c = p(img, w, h, ix  , iy+1);
	float d = p(img, w, h, ix+1, iy+1);
	float r = evaluate_bilinear_cell(a, b, c, d, x-ix, y-iy);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float x = (i - 0.5)/2;
		float y = (j - 0.5)/2;
		out[ow*j+i] = bilinear_interpolation(in, iw, ih, x, y);
	}
}

// THE FILL_BILL ALGORITHM
// @y: output image
// @x: input image, whose NAN pixels are to be filled
// @w: image width
// @h: image height
// @n: number of recursive steps
void fill_bill_recursive(float *y, float *x, int w, int h, int n)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (n > 1)
	{
		// the "s" suffix means "small"
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *xs = xmalloc(ws * hs * sizeof*xs);
		float *ys = xmalloc(ws * hs * sizeof*ys);
		zoom_out_by_factor_two(xs, ws, hs, x, w, h);
		fill_bill_recursive(ys, xs, ws, hs, n - 1);
		zoom_in_by_factor_two(init, w, h, ys, ws, hs);
		free(xs);
		free(ys);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	for (int i = 0; i < w*h; i++)
		y[i] = isfinite(x[i]) ? x[i] : init[i];
	free(init);
}

// run fill_bill at each channel of a multi-channel image
void fill_bill_split(float *out, float *in, int w, int h, int pd)
{
	int max_scales = 100; // 100 scales ought to be enough for anybody
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		fill_bill_recursive(outl, inl, w, h, max_scales);
	}
}


// ddp the following two lines to enable/disable the "main" function
#undef MAIN_FILL_BILL
#define MAIN_FILL_BILL

#ifdef MAIN_FILL_BILL
#include "iio.h"
int main(int argc, char *argv[])
{
	if (argc != 4)
		return fprintf(stderr, "usage:\n\t"
				"%s data.png mask.png out.png\n", *argv);
				//0 1        2        3
	char *filename_in = argv[1];
	char *filename_mask = argv[2];
	char *filename_out = argv[3];

	int w[2], h[2], pd;
	float *in = iio_read_image_float_split(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	float *out = xmalloc(*w**h*pd*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			for (int l = 0; l < pd; l++)
				in[*w**h*l+i] = NAN;

	fill_bill_split(out, in, *w, *h, pd);

	iio_save_image_float_split(filename_out, out, *w, *h, pd);

	return 0;
}
#endif//MAIN_FILL_BILL
