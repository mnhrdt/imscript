#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static bool global_fold = 0;

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

void ihough(float *transform, int n_theta, int n_rho, float *imag, int w, int h,
		float minmag)
{
	// fill counts to zero
	for (int i = 0; i < n_theta * n_rho; i++)
		transform[i] = 0;

	// for each bubble
	// accumulate the corresponding Hough curve by its magnitude
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float x = i - w/2;
		float y = j - h/2;
		float magnitude = imag[j*w+i];

		if (magnitude < minmag)
			continue;

		for (int i_theta = 0; i_theta < n_theta; i_theta++)
		{
			float ti_theta = i_theta;
			float theta = (2.0 * M_PI * ti_theta) / n_theta;
			float alpha = M_PI + atan2(y, x);
			float rho = hypot(x, y) * cos(alpha - theta);
			if (global_fold && ti_theta > n_theta/2)
			{
				ti_theta -= n_theta/2.0;
				rho *= -1;
			}
			//if (rho < 0) {
			//	rho = -rho;
			//	ti_theta = ti_theta - n_theta/2;
			//}
			//if (ti_theta > n_theta/2)
			//	ti_theta = ti_theta - n_theta/2;
			int i_rho = n_rho   * (0.5 + rho   / hypot(w, h) );
			//if (global_fold) ti_theta *= 2;
			if (insideP(n_theta, n_rho, (int)ti_theta, i_rho))
				transform[i_rho*n_theta + (int)ti_theta] += magnitude;
		}
	}
}


#define MAIN_IHOUGH2

#ifdef MAIN_IHOUGH2
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	global_fold = pick_option(&c, &v, "f", NULL);
	// process input arguments
	if (c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s nrho ntheta mtres <intensity >hough\n", *v);
		//               0  1    2      3
		return 1;
	}
	int nrho = atoi(v[1]);
	int ntheta = atoi(v[2]);
	double mtres = atof(v[3]);

	// read input gradient
	int w, h, pd;
	float *bubbles = iio_read_image_float_vec("-", &w, &h, &pd);
	if (pd != 1) return fprintf(stderr, "I expect an intensity!\n");

	// compute transform
	float *transform = malloc(nrho * ntheta * sizeof*transform);
	ihough(transform, nrho, ntheta, bubbles, w, h, mtres);

	// save output image
	//if (global_fold) ntheta /= 2;
	iio_save_image_float_vec("-", transform, nrho, ntheta, 1);

	// cleanup and exti
	free(bubbles);
	free(transform);
	return 0;
}
#endif//MAIN_IHOUGH2
