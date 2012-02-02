// produce a correlated sequence of random vector fields
// SIGMA = gaussian strength of field
// ETA = gaussian strength of dependence
// PHI = gaussian strenght of the temporal correlation
//
// TODO: replace ETA and PHI by a more general 3d matrix

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

//void fill_random_field(float *f, int w, int h, float sigma, float eta)
//{
//	for (int i = 0; i < w * h * 2; i++)
//		f[i] = sigma*random_normal();
//
//	void gblur(float *y, float *x, int w, int h, int pd, float s);
//	gblur(f, f, w, h, 2, eta);
//
//}


void fill_random_fields(float *f, int w, int h, int d,
					float sigma, float eta, float tau)
{
	for (int i = 0; i < w * h * d * 2; i++)
		f[i] = sigma*random_normal();

	float s[3] = {eta, eta, tau};
	void gblur3d(float *,float *,int, int, int, int, float s[3]);
	gblur3d(f, f, w, h, d, 2, s);
}

int main(int c, char *v[])
{
	if (c != 8) {
		fprintf(stderr, "usage:\n\t%s w h d sigma eta phi opat\n", *v);
		//                          0 1 2 3 4     5   6   7
		return EXIT_FAILURE;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int d = atoi(v[3]);
	float sigma = atof(v[4]);
	float eta = atof(v[5]);
	float phi = atof(v[6]);
	char *out_pattern = v[7];

	float *f = xmalloc(w * h * d * 2 * sizeof*f);

	fill_random_fields(f, w, h, d, sigma, eta, phi);

	for (int i = 0; i < d; i++) {
		char buf[0x200];
		snprintf(buf, 0x200, out_pattern, i);
		iio_save_image_float_vec(buf, f+(w*h*2)*i, w, h, 2);
	}

	free(f);

	return EXIT_SUCCESS;
}
