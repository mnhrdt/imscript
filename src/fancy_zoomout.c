#include <math.h>
#include <assert.h>
#include "fancy_image.h"

static double average_of_non_nans(double *x, int n)
{
	int cx = 0;
	double ax = 0;
	for (int i = 0; i < n; i++)
	{
		if (isnan(x[i])) continue;
		cx += 1;
		ax += x[i];
	}
	return cx ? ax / cx : NAN;
}

static double combine_4doubles(double v[4], int m)
{
	if (m == 'f') return v[0];
	if (m == 'v') return (v[0]+v[1]+v[2]+v[3])/4;
	if (m == 'V') return average_of_non_nans(v, 4);
	if (m == 'i') return fmin(fmin(v[0],v[1]), fmin(v[2],v[3]));
	if (m == 'a') return fmax(fmax(v[0],v[1]), fmax(v[2],v[3]));
	return m;//NAN;
}

static void fancy_zoom_out_by_factor_two(char *fname_out, char *fname_in, int m)
{
	// open input image
	struct fancy_image *a = fancy_image_open(fname_in, "r");

	// read information from input image
	int tw = 0, th = 0, fmt = 0, bps = 0;
	int tiffo = fancy_image_leak_tiff_info(&tw, &th, &fmt, &bps, a);

	// create output image of the appropriate size and options
	int pw = ceil(a->w / 2.0);
	int ph = ceil(a->h / 2.0);
	struct fancy_image *b = fancy_image_create(fname_out,
			"w=%d,h=%d,pd=%d,bps=%d,fmt=%d,tw=%d,th=%d",
			pw, ph, a->pd, bps, fmt, tw, th);

	// fill-in the zoomed-out image
	for (int j = 0; j < b->h; j++)
	for (int i = 0; i < b->w; i++)
	for (int l = 0; l < b->pd; l++)
	{
		int ii = 2 * i;
		int jj = 2 * j;
		double v[4] = {
			fancy_image_getsample(a, ii + 0, jj + 0, l),
			fancy_image_getsample(a, ii + 1, jj + 0, l),
			fancy_image_getsample(a, ii + 0, jj + 1, l),
			fancy_image_getsample(a, ii + 1, jj + 1, l)
		};
		double r = combine_4doubles(v, m);
		fancy_image_setsample(b, i, j, l, r);
	}

	// close both images
	fancy_image_close(b);
	fancy_image_close(a);
}

#include <stdio.h>
int main(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s {f|v|i|a} in.tiff out.tiff\n", *v);
		//               0   1        2       3
		return 1;
	}
	int op = v[1][0];
	char *filename_in = v[2];
	char *filename_out = v[3];

	fancy_zoom_out_by_factor_two(filename_out, filename_in, op);

	return 0;
}
