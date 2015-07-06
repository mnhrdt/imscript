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

void ghough(float *transform, int n_theta, int n_rho, float *grad, int w, int h)
{
	// fill counts to zero
	for (int i = 0; i < n_theta * n_rho; i++)
		transform[i] = 0;

	// for each gradient vector,
	// increase the corresponding Hough count by its norm
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float *g = grad + 2 * (j*w + i);
		float x = i - w/2;
		float y = j - h/2;
		float ng = hypot(g[0], g[1]);
		float theta = M_PI + atan2(g[1], g[0]);
		float alpha = M_PI + atan2(y, x);
		float rho = hypot(x, y) * cos(alpha - theta);
		float ti_theta = n_theta * theta / (2 * M_PI);
		if (global_fold && ti_theta > n_theta/2)
		{
			ti_theta -= n_theta/2.0;
			rho *= -1;
		}
		int i_rho   = n_rho   * (0.5 + rho   / hypot(w, h) );
		if (insideP(n_theta, n_rho, (int)ti_theta, i_rho))
			transform[i_rho*n_theta + (int)ti_theta] += ng;
	}
}


#define MAIN_GHOUGH2

#ifdef MAIN_GHOUGH2
#include <stdio.h>
#include <stdlib.h>
#include "pickopt.c"
#include "iio.h"
int main(int c, char *v[])
{
	global_fold = pick_option(&c, &v, "f", NULL);
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
				"%s nrho ntheta <gradient >hough\n", *v);
		//               0  1    2
		return 1;
	}
	int nrho = atoi(v[1]);
	int ntheta = atoi(v[2]);

	// read input gradient
	int w, h, pd;
	float *gradient = iio_read_image_float_vec("-", &w, &h, &pd);
	if (pd != 2) return fprintf(stderr, "I expect a gradient!\n");

	// compute transform
	float *transform = malloc(nrho * ntheta * sizeof*transform);
	ghough(transform, nrho, ntheta, gradient, w, h);

	// save output image
	iio_save_image_float_vec("-", transform, nrho, ntheta, 1);

	// cleanup and exti
	free(gradient);
	free(transform);
	return 0;
}
#endif//MAIN_GHOUGH2
