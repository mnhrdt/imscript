// input: a complex-valued image
// output: its square root, unwrapped by ising

#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include "random.c"
#include "iio.h"

// distance of two angles on the unit circle
static float angular_distance(float a, float b)
{
	a -= b;
	while( a <= -M_PI ) a += 2*M_PI;
	while( a >   M_PI ) a -= 2*M_PI;
	return fabs(a) / M_PI;
}

static
float local_energy(complex float *y, int w, int h, int i, int j, bool flip)
{
	// extract values of neighbors
	int nn = 0;
	complex float n[4];
	if (i > 0)   n[nn++] = y[ j   *w + i-1];
	if (j > 0)   n[nn++] = y[(j-1)*w + i  ];
	//if (i < w-1) n[nn++] = y[ j   *w + i+1];
	//if (j < h-1) n[nn++] = y[(j+1)*w + i  ];

	// extract value at the current site
	complex float b = y[j*w+i] * (flip ? -1 : 1);

	// compute laplacian energy energy
	float r = 0;
	for (int k = 0; k < nn; k++)
		r += cabs(b - n[k]);//angular_distance(carg(b), carg(n[k]));
	return r;
}

static void sqroot_ising(complex float *y, complex float *x, int w, int h,
		float temperature, float nsteps)
{
	// 1. compute an arbitrarily chosen square root
	for (int i = 0; i < w*h; i++)
		y[i] = csqrt(x[i]);

	// 2. initialize with random signs
	for (int i = 0; i < w*h; i++)
		if (random_uniform() > 0.5)
			y[i] *= -1;

	// debugging stuff
	iio_write_image_float_vec("/tmp/ising_0.tiff", (void*)y, w, h, 2);
	float *tmp_1 = malloc(w*h*sizeof(float));
	float *tmp_2 = malloc(w*h*sizeof(float));
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		tmp_1[j*w+i] = local_energy(y, w, h, i, j, true);
		tmp_2[j*w+i] = local_energy(y, w, h, i, j, false);
	}
	iio_write_image_float_vec("/tmp/E1.tiff", tmp_1, w, h, 1);
	iio_write_image_float_vec("/tmp/E2.tiff", tmp_2, w, h, 1);

	// 3. metropolis algorithm
	float beta = 1 / temperature;
	for (int k = 0; k < nsteps*w*h; k++)
	{
		//int i = randombounds(0, w-1);
		//int j = randombounds(0, h-1);
		int i = k % w;
		int j = (k / w)%h;
		//if (k%2) i = w - i - 1;
		//if (k%2) j = h - j - 1;
		float E1 = local_energy(y, w, h, i, j, true);
		float E2 = local_energy(y, w, h, i, j, false);
		float dE = E1 - E2;
		if (dE < 0 || random_uniform() <= exp(-beta * dE)) 
			y[j*w+i] *= -1;
	}
}

int main(int c, char *v[])
{
	if (c < 3)
		fprintf(stderr, "usage:\n\t%s temp nsteps [in [out]]\n", *v);
	//                                  0 1    2       3   4
	float temperature = atof(v[1]);
	float nsteps = atof(v[2]);
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (pd != 2)
		return fprintf(stderr, "ERROR: "
				"input should be complex (got %d)\n", pd);
	complex float *y = malloc(w*h*sizeof*y);
	sqroot_ising(y, (void*)x, w, h, temperature, nsteps);
	iio_write_image_float_vec(filename_out, (void*)y, w, h, 2);
	free(x);
	free(y);
	return 0;
}
