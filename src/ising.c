// simplest ising model (with metropolis sampling)

#include <assert.h>
#include "random.c"

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

// s: image of signs (initialized, updated by thus function)
// J: coupling image
// H: exterior field (optional, to give a background context)
// temperature: strictly positive
// nsteps: number of steps to sample
void ising_metropolis(float *s, float *J, float *H, int w, int h,
		float temperature, float nsteps)
{
	// assert input consistency
	for (int i = 0; i < w*h; i++)
		assert(abs(abs(s[i])-1) < 1e-9);
	assert(temperature > 0);
	assert(nsteps > 0);

	// metropolis algorithm
	double beta = 1 / temperature;
	for (int k = 0; k < nsteps; k++)
	{
		int i = randombounds(0, w - 1);
		int j = randombounds(0, h - 1);
		float bn = local_field(s, J, H, w, h, i, j);
		float P1bn = 1 / (1 + exp(-2 * beta * bn)); // P(s = +1 | bn)
		float U = random_uniform();
		s[j*w+i] = U < P1bn ? 1 : -1;
	}
}
