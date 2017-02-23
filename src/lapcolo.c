// laplace-beltrami colorizer

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OMIT_LABPEDIAG_MAIN
#include "lapbediag.c"


#include "smapa.h"
SMART_PARAMETER(TSTEP,0.1)
SMART_PARAMETER(NITER,1)
SMART_PARAMETER(NSCALES,1)

static bool grayP(float *g, int n)
{
	for (int i = 1; i < n; i++)
		if (g[i] != g[0])
			return false;
	return true;
}

static float norm(float *x, int n)
{
	return n ? hypot(*x, norm(x+1,n-1)) : 0;
}

// y gets the hue of x and the intensity of c
static void normalize_intensity(float *y, float *x, float *c, int n)
{
	float mx = norm(x,n);
	float mc = norm(c,n);
	//float mx = 0, mc = 0;
	//for (int i = 0; i < n; i++)
	//{
	//	mc += c[i]/n;
	//	mx += x[i]/n;
	//}
	for (int i = 0; i < n; i++)
		y[i] = mx ? (mc/mx)*x[i] : c[i];
}

void lapbe_colorizer(float *outhue, float *outint,
		float *metric, float *color, int w, int h, int pd)
{
	float tstep = TSTEP();
	int niter = NITER();
	int nscales = NSCALES();

	float *tmpi = xmalloc(w*h*sizeof*tmpi);
	float *tmpo = xmalloc(w*h*sizeof*tmpo);
	for (int l = 0; l < pd; l++)
	{
		for (int i = 0; i < w*h; i++)
			tmpi[i] = grayP(color + i*pd, pd) ? NAN : color[i*pd+l];
		lapbediag_rec(tmpo, metric, tmpi, w, h, tstep, niter, nscales);
		for (int i = 0; i < w*h; i++)
			outhue[i*pd+l] = tmpo[i];
	}
	free(tmpo);
	free(tmpi);

	for (int i = 0; i < w*h; i++)
		normalize_intensity(outint+i*pd, outhue+i*pd, color+i*pd, pd);
}


#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 5) {
		fprintf(stderr, "usage:\n\t"
		"%s metric colors outk outi\n", *argv);
		//0 1      2      3    4
		return 1;
	}
	char *filename_metric = argv[1];
	char *filename_colors = argv[2];
	char *filename_outhue = argv[3];
	char *filename_outint = argv[4];

	int w[2], h[2], pd;
	float *metric = iio_read_image_float(filename_metric, w, h);
	float *colors = iio_read_image_float_vec(filename_colors, w+1, h+1,&pd);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "input image files sizes mismatch");
	float *outhue = xmalloc(pd**w**h*sizeof*outhue);
	float *outint = xmalloc(pd**w**h*sizeof*outhue);


	lapbe_colorizer(outhue, outint, metric, colors, *w, *h, pd);

	iio_write_image_float_vec(filename_outhue, outhue, *w, *h, pd);
	iio_write_image_float_vec(filename_outint, outint, *w, *h, pd);

	free(outhue); free(outint); free(metric); free(colors);
	return 0;
}
