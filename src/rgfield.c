// produce a random vector field
// SIGMA = gaussian strength of field
// ETA = gaussian strength of dependence

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"

#include "xmalloc.c"

#include "smapa.h"
SMART_PARAMETER_SILENT(RSEED,0)

static double random_uniform(void)
{
	static bool firstrun = true;
	if (firstrun) {
		firstrun = false;
		srand(RSEED());
	}
	return rand()/(RAND_MAX+1.0);
}

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y1 = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y1;
}

void fill_random_field(float *f, int w, int h, float sigma, float eta)
{
	for (int i = 0; i < w * h * 2; i++)
		f[i] = sigma*random_normal();

	void gblur(float *y, float *x, int w, int h, int pd, float s);
	gblur(f, f, w, h, 2, eta);

}


//void fill_random_fields(float *f, int w, int h, int t,
//					float sigma, float eta, float tau)
//{
//	for (int i = 0; i < w * h * t * 2; i++)
//		f[i] = sigma*random_normal();
//
//	gblur(f, f, w, h, 2, eta);
//}

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s w h sigma eta\n", *v);
		//                          0 1 2 3     4
		return EXIT_FAILURE;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	float sigma = atof(v[3]);
	float eta = atof(v[4]);

	float *f = xmalloc(w * h * 2 * sizeof*f);

	fill_random_field(f, w, h, sigma, eta);

	iio_save_image_float_vec("-", f, w, h, 2);

	free(f);

	return EXIT_SUCCESS;
}
