// local lgblur (a section of the scale space)
// this version is optimized por generating multiple versions of the same image
// usage: lgblur2 infile.png varpat outpat FIRST LAST


#include <assert.h>
#include <math.h>
#include <stdio.h>


static float index_to_scale(int n, float s0, float sf, int i)
{
	float r = -1;
	if (i < 0) r = s0;
	else if (i >= n) r = sf;
	//else r = s0 * pow(sf/s0, i/(n-1.0));
	else r = s0 + i * (sf - s0)/(n - 1);
	assert(r >= s0);
	assert(r <= sf);
	return r;
}

static int scale_to_index(int n, float s0, float sf, float s)
{
	int r = -1;
	if (s <= s0) r = 0;
	else if (s >= sf) r = n - 1;
	//else r = lround((n - 1)*log(s/s0)/log(sf/s0));
	else r = lround((n - 1)*(s-s0)/(sf-s0));
	assert(r >= 0);
	assert(r < n);
	return r;
}

#define OMIT_GBLUR_MAIN
#include "gblur.c"

// utility headers used only in the "main" function
#include <stdio.h>

#include "iio.h"
#include "fail.c"
#include "xmalloc.c"

// fill a pyramid of images of the same size
void build_gaussian_pyramid(float *p, int npyr, float sfirst, float slast,
		float *x, int w, int h, int pd)
{
	assert(0 < sfirst);
	assert(sfirst < slast);
	for (int i = 0; i < npyr; i++)
	{
		float s = index_to_scale(npyr, sfirst, slast, i);
		fprintf(stderr, "s[%d] = %g\n", i, s);
		gblur(p+i*w*h*pd, x, w, h, pd, s);

		char buf[FILENAME_MAX];
		snprintf(buf, FILENAME_MAX, "/tmp/pyra_%03d.png", i);
		iio_save_image_float_vec(buf, p+i*w*h*pd, w, h, pd);
	}
}

void apply_local_blur(float *oy, float *sigma, float *px,
		int npyr, float sfirst, float slast, int w, int h, int pd)
{
	float (*p)[h][w][pd] = (void*)px;
	float (*y)[w][pd] = (void*)oy;
	for (int l = 0; l < pd; l++)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float s = sigma[j*w+i];
		int sidx = scale_to_index(npyr, sfirst, slast, s);
		y[j][i][l] = p[sidx][j][i][l];
	}
}

#include "bicubic.c"

void apply_general_model(float *oy, float *sigma, float *uv, float *px,
		int npyr, float sfirst, float slast, int w, int h, int pd)
{
	float (*p)[h][w][pd] = (void*)px;
	float (*y)[w][pd] = (void*)oy;
	float (*f)[w][2] = (void*)uv;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float s = sigma[j*w+i];
		int sidx = scale_to_index(npyr, sfirst, slast, s);
		float dx = uv[j][i][0];
		float dy = uv[j][i][1];
		float v[pd];
		bicubic_interpolation(v, p[sidx][0][0], w, h, pd, i+dx, j+dy);
		for (int l = 0; l < pd; l++)
			y[j][i][l] = v[l];
	}
	for (int l = 0; l < pd; l++)

	{
	}
}

int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t"
			"%s in.png varpat uvpat outpat FIRST LAST\n", *v);
		//        0 1      2      3     4      5     6
		return 1;
	}
	char *filename_in = v[1];
	char *filepattern_var = v[2];
	char *filepattern_flo = v[3];
	char *filepattern_out = v[4];
	int first = atoi(v[5]);
	int last = atoi(v[6]);

	int npyr = 60;
	float sfirst = 0.0001;
	float slast = 30;

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *px = xmalloc(w*h*pd*npyr*sizeof*px);
	build_gaussian_pyramid(px, npyr, sfirst, slast, x, w, h, pd);
	float *y = xmalloc(w*h*pd*sizeof*y);

	for (int i = first; i < last; i++)
	{
		char filename_var[FILENAME_MAX];
		char filename_flo[FILENAME_MAX];
		char filename_out[FILENAME_MAX];
		snprintf(filename_var, FILENAME_MAX, filepattern_var, i);
		snprintf(filename_flo, FILENAME_MAX, filepattern_flo, i);
		snprintf(filename_out, FILENAME_MAX, filepattern_out, i);

		int ww, hh, ppdd;
		float *sigma = iio_read_image_float(filename_var, &ww, &hh);
		if (w != ww || h != hh)
			fail("variances and image sizes mismatch");
		float *uv = iio_read_image_float(filename_flo, &ww, &hh, &ppdd);
		if (ppdd != 2)
			fail("expect vector field on file \"%s\"",filename_flo);
		if (w != ww || h != hh)
			fail("flow and image sizes mismatch");

		apply_general_model(y, sigma,uv, px, npyr,sfirst,slast, w,h,pd);
		//apply_local_blur(y, sigma, px, npyr, sfirst, slast, w, h, pd);
		iio_save_image_float_vec(filename_out, y, w, h, pd);
		free(sigma);
	}

	free(x); free(y);
	return 0;
}
