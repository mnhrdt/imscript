#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fail.c"
#include "xmalloc.c"
#include "random.c"
#include "bilinear_interpolation.c"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#include "smapa.h"
SMART_PARAMETER(LIC_IFACTOR,10)
SMART_PARAMETER(LIC_TRES,1)

void line_integral_convolution(float *view, float *flow, int w, int h)
{
	float (*f)[w][2] = (void*)flow;
	float (*v)[w] = (void*)view;

	FORJ(h) FORI(w)
		v[j][i] = random_uniform();

	int niter = LIC_IFACTOR() * w * h;
	float litres = LIC_TRES();
	FORL(niter) {
		int i = randombounds(0, w-1);
		int j = randombounds(0, h-1);
		float fn[2] = {f[j][i][0], f[j][i][1]};
		float fnn = hypot(fn[0], fn[1]);
		if (fnn > litres) {
			fn[0] /= litres*fnn;
			fn[1] /= litres*fnn;
		}
		float fnext[2] = {i+fn[0], j+fn[1]};
		float fprev[2] = {i-fn[0], j-fn[1]};
		float vals[3] = {
			bilinear_interpolation_at(view,w,h, fnext[0],fnext[1]),
			v[j][i],
			bilinear_interpolation_at(view,w,h, fprev[0],fprev[1])};
		float valsum = (vals[0] + 2*vals[1] + vals[2])/4;
		v[j][i] = valsum;
	}
}


#ifndef OMIT_MAIN
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3 && c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [flow [view]]\n", *v);
				//          0  1     2
		return EXIT_FAILURE;
	}
	char *infile = c > 1 ? v[1] : "-";
	char *outfile = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *flow = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) fail("input is not a vector field");

	float *view = xmalloc(w * h * sizeof*view);
	line_integral_convolution(view, flow, w, h);
	iio_save_image_float_vec(outfile, view, w, h, 1);

	free(view);
	free(flow);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
