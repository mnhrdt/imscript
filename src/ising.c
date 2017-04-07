// simplest ising model (with metropolis sampling)

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "random.c"

static void flip_spin (float *x, int w, int h, int i, int j) { x[j*w+i] *= -1; }
static void flip_plus (float *x, int w, int h, int i, int j) { x[j*w+i]  =  1; }
static void flip_minus(float *x, int w, int h, int i, int j) { x[j*w+i]  = -1; }

// NOTE: the commented part below is specific to signed fields
// a slightly more general implementation is given at the end
//
//static float local_energy(float *x, int w, int h, int i, int j, bool flip)
//{
//	// extract values of neighbors
//	int nn = 0;
//	float n[4];
//	if (i > 0)   n[nn++] = y[ j   *w + i-1];
//	if (j > 0)   n[nn++] = y[(j-1)*w + i  ];
//	if (i < w-1) n[nn++] = y[ j   *w + i+1];
//	if (j < h-1) n[nn++] = y[(j+1)*w + i  ];
//
//	// extract value at the current site
//	float b = x[j*w+i] * flip ? -1 : 1;
//
//	// compute laplacian energy
//	float r = 0;
//	for (int k = 0; k < nn; k++)
//		r += b - n[k];
//	return fabs(r);
//}
//
//// computation using the Metropolis Algorithm
//// (alternatives: rejection sampling, Gibbs, annealing)
//void run_ising_sign_step(float *x, int w, int h, float temperature, int nsteps)
//{
//	if (nsteps < w*h)
//		fprintf(stderr, "WARNING: Ising with few steps!\n");
//
//	float beta = 1 / temperature;
//
//	for (int k = 0; k < nsteps; k++)
//	{
//		int i = randombounds(0, w - 1);
//		int j = randombounds(0, h - 1);
//		float E1 = local_energy(x, w, h, i, j, true);
//		float E2 = local_energy(x, w, h, i, j, false);
//		float dE = E1 - E2;
//
//		if (dE < 0)
//			flip_spin(x, w, h, i, j);
//		else
//			if (random_uniform() <= exp(-beta * dE))
//				flip_spin(x, w, h, i, j);
//	}
//}

static float local_field(float *s, float *J, float *H, int w, int h,
		int i, int j)
{
	// extract signs of neighbors
	int nn = 0;
	float n[4];
	if (i > 0)   n[nn++] = s[ j   *w + i-1];
	if (j > 0)   n[nn++] = s[(j-1)*w + i  ];
	if (i < w-1) n[nn++] = s[ j   *w + i+1];
	if (j < h-1) n[nn++] = s[(j+1)*w + i  ];

	// initialize energy counter with exterior field (or with zero)
	float r = 0;
	if (H)
		r += H[j*w+i];

	// accumulate energy
	for (int k = 0; k < nn; k++)
		r += n[k] * J[j*w+i];
	return r;
}

// s: image of signs (initialized, updated by function)
// J: coupling image
// H: exterior field (optional)
// temperature: strictly positive
// nsteps: number of steps to sample
void ising_general(float *s, float *J, float *H, int w, int h,
		float temperature, float nsteps)
{
	// assert input consistency
	for (int i = 0; i < w*h; i++)
		assert(abs(abs(s[i])-1) < 1e-9);
	assert(temperature > 0);
	assert(nsteps > 0);
	if (nsteps < 17*w*h)
		fprintf(stderr, "WARNING: Ising with few steps!\n");

	// parameter
	double beta = 1 / temperature;

	// metropolis iterations
	for (int k = 0; k < nsteps; k++)
	{
		int i = randombounds(0, w - 1);
		int j = randombounds(0, h - 1);
		float bn = local_field(s, J, H, w, h, i, j);
		float P1bn = 1 / (1 + exp(-2 * beta * bn)); // P(s = +1 | bn)
		assert(P1bn >= 0);
		assert(P1bn <= 1);

		float U = random_uniform();
		if (U < P1bn)
			flip_plus(s, w, h, i, j);
		else
		       flip_minus(s, w, h, i, j);
	}
}
