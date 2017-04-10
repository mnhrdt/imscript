// input: a complex-valued image
// output: its square root, unwrapped by ising

#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdbool.h>
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
	complex float n[8];
	//for (int dy = -1; dy <= 1; dy++)
	//for (int dx = -1; dx <= 1; dx++)
	//{
	//	int ii = i + dx;
	//	int jj = j + dy;
	//	if (ii >= 0 && ii < w && jj >= 0 && jj < h)
	//		n[nn++] = y[jj*w+ii];
	//}
	if (i > 0)   n[nn++] = y[ j   *w + i-1];
	//if (i < w-1) n[nn++] = y[ j   *w + i+1];
	if (j > 0)   n[nn++] = y[(j-1)*w + i  ];
	//if (j < h-1) n[nn++] = y[(j+1)*w + i  ];

	// extract value at the current site
	complex float b = y[j*w+i] * (flip ? -1 : 1);

	// compute laplacian energy
	float r = 0;
	for (int k = 0; k < nn; k++)
		//r += cabs(b - n[k]);
		r += angular_distance(carg(b), carg(n[k]));
	return r;
}

// compute the complex square root by Ising unwarping
static void sqroot_ising(complex float *y, complex float *x, int w, int h,
		float temperature, float nsteps, complex float *init)
{
	// 1. compute an arbitrarily chosen square root
	for (int i = 0; i < w*h; i++)
		y[i] = csqrt(x[i]);

	// 2. initialize with random signs (or with "init" signs, if given)
	if (init)
		for (int i = 0; i < w*h; i++)
			if (creal(init[i] * conj(y[i])) < 0)
				y[i] *= -1;
	if (!init)
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

// extrapolate by nearest value (useful for Neumann boundary conditions)
static complex float getpixel_1(complex float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// zoom out by factor two
static void zoom_out_by_factor_two(complex float *out, int ow, int oh,
		complex float *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		complex float a[4], m = 0;
		a[0] = getpixel_1(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel_1(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel_1(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel_1(in, iw, ih, 2*i+1, 2*j+1);
		for (int k = 0; k < 4; k++)
			m += a[k];
		out[ow*j + i] = m;// / 4;
	}
}

// evaluate a bilinear cell at the given point
static complex float evaluate_bilinear_cell(complex float a[4],
		float x, float y)
{
	complex float r = 0;
	r += a[0] * (1-x) * (1-y);
	r += a[1] * ( x ) * (1-y);
	r += a[2] * (1-x) * ( y );
	r += a[3] * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static complex float bilinear_interpolation(complex float *x, int w, int h,
		float p, float q)
{
	int ip = floor(p); // note: a cast to int fails when p<0
	int iq = floor(q);
	complex float a[4];
	a[0] = getpixel_1(x, w, h, ip  , iq  );
	a[1] = getpixel_1(x, w, h, ip+1, iq  );
	a[2] = getpixel_1(x, w, h, ip  , iq+1);
	a[3] = getpixel_1(x, w, h, ip+1, iq+1);
	complex float r = evaluate_bilinear_cell(a, p-ip, q-iq);
	return r;
}


// zoom-in by replicating pixels into 2x2 blocks
static void zoom_in_by_factor_two(complex float *out, int ow, int oh,
		complex float *in, int iw, int ih)
{
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		complex float x = (i - 0.5) / 2;
		complex float y = (j - 0.5) / 2;
		out[ow*j+i] = bilinear_interpolation(in, iw, ih, x, y);
		//out[ow*j+i] = getpixel_1(in, iw, ih, i/2, j/2);
	}
}



static void pyrsquising_rec(complex float *y, complex float *x, int w, int h)
{
	// compute initialization by multi-scale
	complex float *init = malloc(w*h*sizeof*init);
	if (w > 1 || h > 1) {
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		complex float *xs = malloc(ws * hs * sizeof*xs);
		complex float *ys = malloc(ws * hs * sizeof*ys);
		zoom_out_by_factor_two(xs, ws, hs, x, w, h);
		pyrsquising_rec(ys, xs, ws, hs);
		zoom_in_by_factor_two(init, w, h, ys, ws, hs);
		free(xs);
		free(ys);
	} else { // size 1x1
		*init = csqrt(*x);
	}

	// refine the low-res initialization
	sqroot_ising(y, x, w, h, 3.0/sqrt(w*h), 3, init);
	//sqroot_ising(y, x, w, h, 0.023, 243, init);

	free(init);
}

// pyramidally unwarped square root
static void pyramidal_sqising(complex float *y, complex float *x, int w, int h)
{
	pyrsquising_rec(y, x, w, h);
}


#include "pickopt.c"
int main(int c, char *v[])
{
	char *filename_init = pick_option(&c, &v, "i", "");
	if (c < 3)
		return fprintf(stderr, "usage:\n\t"
				"%s [-i init] temp nsteps [in [out]]\n", *v);
	//                        0           1    2       3   4
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
	float *init = NULL;
	if (*filename_init) {
		int ww, hh;
		init = iio_read_image_float_vec(filename_init, &ww, &hh, &pd);
		if (ww != w || hh != h || pd != 2)
			return fprintf(stderr, "ERORR: bat init\n");
		fprintf(stderr, "using init from file \"%s\"\n", filename_init);
	}
	if (temperature > 0)
		sqroot_ising(y,(void*)x, w,h, temperature,nsteps, (void*)init);
	else
		pyramidal_sqising(y, (void*)x, w, h);
	iio_write_image_float_vec(filename_out, (void*)y, w, h, 2);
	free(x);
	free(y);
	if (init) free(init);
	return 0;
}
