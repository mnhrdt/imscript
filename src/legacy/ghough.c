#include <assert.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1416
#endif

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
		float ng = hypot(g[0], g[1]);
		float theta = M_PI + atan2(g[1], g[0]);
		float alpha = M_PI + atan2(j, i);
		float rho = hypot(w,h) + hypot(i, j) * cos(alpha - theta);
		float max_rho = 2*hypot(w, h);
		float max_theta = 2*M_PI;
		int i_rho   = n_rho   * rho   / max_rho;
		int i_theta = n_theta * theta / max_theta;
		assert(0 <= i_rho);
		assert(0 <= i_theta);
		assert(i_rho < n_rho);
		assert(i_theta < n_theta);
		transform[i_rho*n_theta + i_theta] += ng;
	}
}


#define MAIN_GHOUGH

#ifdef MAIN_GHOUGH
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
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
	iio_write_image_float_vec("-", transform, nrho, ntheta, 1);

	// cleanup and exti
	free(gradient);
	free(transform);
	return 0;
}
#endif//MAIN_GHOUGH
