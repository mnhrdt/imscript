// Produce a correlated sequence of random vector fields
//
// SIGMA = gaussian strength of field
//
// Entries of spatiotemporal (x,y,t) covariance matrix:
//
// E1 E2 E3
//    E4 E5
//       E6
//

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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

#define OMIT_GBLUR_MAIN
#include "gblur.c"

void fill_random_fields(float *f, int w, int h, int d,
		float sigma, float eta[6])
{
	for (int i = 0; i < w * h * d * 2; i++)
		f[i] = sigma*random_normal();

	gblur3dm(f, f, w, h, d, 2, eta);
}

int main(int c, char *v[])
{
	if (c != 12) {
		fprintf(stderr, "usage:\n\t%s w h d sigma "
		//                          0 1 2 3 4
				"e1 e2 e3 e4 e5 e6 opat\n", *v);
		//               5  6  7  8  9  10 11
		return EXIT_FAILURE;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int d = atoi(v[3]);
	float sigma = atof(v[4]);
	float eta[6];
	for (int i = 0; i < 6; i++)
		eta[i] = atof(v[5+i]);
	char *out_pattern = v[11];

	float *f = xmalloc(w * h * d * 2 * sizeof*f);

	fill_random_fields(f, w, h, d, sigma, eta);

	for (int i = 0; i < d; i++) {
		char buf[0x200];
		snprintf(buf, 0x200, out_pattern, i);
		iio_write_image_float_vec(buf, f+(w*h*2)*i, w, h, 2);
	}

	free(f);

	return EXIT_SUCCESS;
}
