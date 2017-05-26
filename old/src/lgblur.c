// local lgblur (a section of the scale space)
// usage: lgblur variances [in [out]]


#include <math.h>
#include <stdio.h>

static float unnormalized_gaussian_function(float sigma, float x)
{
	return exp(-x*x/(2*sigma*sigma));
}

#define KWMAX 39

static int gaussian_kernel_width(float sigma)
{
	float radius = 3 * fabs(sigma);
	int r = ceil(1 + 2*radius);
	if (r < 1) r = 1;
	if (r > KWMAX) r = KWMAX;
	return r;
}

static void fill_gaussian_kernel(float *k, int w, float s)
{
	int c = (w - 1)/2;
	for (int j = 0; j < w; j++)
	for (int i = 0; i < w; i++)
		k[j*w+i] = unnormalized_gaussian_function(s, hypot(i-c,j-c));
//
//	float m = 0;
//	for (int i = 0; i < w*w; i++)
//		m += k[i];
//	for (int i = 0; i < w*w; i++)
//		k[i] /= m;
}

static void fill_square(float *s, int ns, float *x, int w, int h, int i, int j)
{
	for (int q = 0; q < ns; q++)
	for (int p = 0; p < ns; p++)
	{
		int ii = i + p - (ns - 1)/2;
		int jj = j + q - (ns - 1)/2;
		float g = NAN;
		if (ii >= 0 && jj >= 0 && ii < w && jj < h)
			g = x[jj*w+ii];
		s[q*ns+p] = g;
	}
}

float local_gaussian_blur_at(float s, float *x, int w, int h, int i, int j)
{
	if (s < 0.1) return x[j*w+i];
	int kw = gaussian_kernel_width(s);
	float k[kw*kw], c[kw*kw];
	fill_gaussian_kernel(k, kw, s);
	fill_square(c, kw, x, w, h, i, j);

	float m = 0;
	for (int l = 0; l < kw*kw; l++)
		if (isfinite(c[l]))
			m += k[l];
	for (int l = 0; l < kw*kw; l++)
		k[l] /= m;

	float r = 0;
	for (int l = 0; l < kw*kw; l++)
		if (isfinite(c[l]))
			r += k[l]*c[l];

	return r;
}

void local_gaussian_blur(float *y, float *s, float *x, int w, int h, int pd)
{
	for (int l = 0; l < pd; l++)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[w*h*l+(j*w+i)]
			= local_gaussian_blur_at(s[j*w+i], x+w*h*l, w,h, i,j);
}


// utility headers used only in the "main" function
#include <stdio.h>

#include "iio.h"
#include "fail.c"
#include "xmalloc.c"

int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s variances [in [out]]\n", *v);
		//                          0 1          2   3
		return 1;
	}
	char *filename_sigma = v[1];
	char *filename_in = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	int w, h, ww, hh, pd;
	float *sigma = iio_read_image_float(filename_sigma, &w, &h);
	float *x = iio_read_image_float_split(filename_in, &ww, &hh, &pd);
	if (w != ww || h != hh) fail("variances and image size mismatch");

	float *y = xmalloc(w*h*pd*sizeof*y);
	local_gaussian_blur(y, sigma, x, w, h, pd);

	iio_write_image_float_split(filename_out, y, w, h, pd);
	free(x); free(y); free(sigma);
	return 0;
}
