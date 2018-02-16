#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

typedef bool (pixel_mask_f)(float*,int,int,int, int,int);

// instance of pixel_mask_f
static
bool pixel_has_nan_sample(float *x, int w, int h, int pd, int i, int j)
{
	for (int l = 0; l < pd; l++)
		if (isnan(x[(j*w+i)*pd+l]))
			return true;
	return false;
}

// instance of pixel_mask_f
static
bool pixel_has_nonpositive_sample(float *x, int w, int h, int pd, int i, int j)
{
	for (int l = 0; l < pd; l++)
	{
		float g = x[(j*w+i)*pd+l];
		if (isnan(g) || g < 0)
			return true;
	}
	return false;
}

#define BAD_MIN(x,y) ((x)<(y)?(x):(y))
#define BAD_MAX(x,y) ((x)>(y)?(x):(y))

void autotrim(float *y, int *out_w, int *out_h, float *x, int w, int h, int pd,
		pixel_mask_f *badP)
{
	// find bounding box
	int i_first = w - 1;
	int j_first = h - 1;
	int i_last  = 0;
	int j_last  = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	if (!badP(x, w, h, pd, i, j))
	{
		i_first = BAD_MIN(i_first, i);
		j_first = BAD_MIN(j_first, j);
		i_last  = BAD_MAX(i_last , i);
		j_last  = BAD_MAX(j_last , j);
	}

	// do the crop
	int ow = i_last - i_first;
	int oh = j_last - j_first;
	fprintf(stderr, "trim if jf il jl ow oh %d %d %d %d %d %d\n",
			i_first, j_first, i_last, j_last, ow, oh);
	assert(ow <= w);
	assert(oh <= h);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < pd; l++)
	{
		int ii = i + i_first;
		int jj = j + j_first;
		y[(j*ow+i)*pd+l] = x[(jj*w+ii)*pd+l];
	}
	*out_w = ow;
	*out_h = oh;
}

#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4)
		return fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
	//                                         0  1   2
	char *filename_in  = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";
	pixel_mask_f *criterion_bad = pixel_has_nonpositive_sample;

	int w, h, pd;//, ow, oh;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *y = malloc(w*h*pd*sizeof*y);
	autotrim(y, &w, &h, x, w, h, pd, criterion_bad);
	iio_write_image_float_vec(filename_out, y, w, h, pd);
	return 0;
}
