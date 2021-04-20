// various algorithms for "content-aware" image-width reducing

#include <assert.h>
#include <math.h>   // hypot, fmin
#include <stdlib.h> // qsort
#include <stdio.h>
#include "iio.h"
#define xmalloc malloc

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

static float sqr(float x)
{
	return x*x;
}

static void compute_energy_field(float *e, float *x, int w, int h, int pd)
{
	float (*X)[w][pd] = (void*)x;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float E = 0;
		if (i==0 || j==0 || i==w-1 || j==h-1)
			E = 0; // protect borders
		else
			for (int l = 0; l < pd; l++)
			{
				// TODO: do something intelligent here
				float x00 = X[j+0][i+0][l];
				float x10 = X[j+1][i+0][l];
				float x01 = X[j+0][i+1][l];
				float xm0 = X[j-1][i+0][l];
				float x0m = X[j+0][i-1][l];
				E += sqr(x10 - x00);
				E += sqr(xm0 - x00);
				E += sqr(x0m - x00);
				E += sqr(x01 - x00);
				//E = fmax(E, fabs(x10 - x00));
				//E = fmax(E, fabs(xm0 - x00));
				//E = fmax(E, fabs(x0m - x00));
				//E = fmax(E, fabs(x01 - x00));
			}
		e[j*w+i] = sqrt(E);
		//e[j*w+i] = E;
	}
}

void reduce_width_naive(float *y, int w2, float *x, int w, int h, int pd)
{
	assert(w2 > 0);
	assert(w2 < w);
	float *e = xmalloc(w * h * sizeof*e); // local energy field
	compute_energy_field(e, x, w, h, pd);
	iio_write_image_float_vec("/tmp/energy_field.npy", e, w, h, 1);

	// accumulate columnar energy
	float c[w];
	for (int i = 0; i < w; i++)
	{
		c[i] = 0;
		for (int j = 0; j < h; j++)
			c[i] += e[j*w+i];
			//c[i] = fmax(c[i], e[j*w+i]);
	}
	//for (int i = 0; i < w; i++)
	//	fprintf(stderr, "%g\n", c[i]/h);

	// sort columns by energy
	float sc[w]; for (int i = 0; i < w; i++) sc[i] = c[i];
	qsort(sc, w, sizeof*c, compare_floats);
	float cutoff = sc[w-w2-1]; // notice: some values may be repeated

	// copy image with lowest energy columns removed
	int ii = 0;
	for (int i = 0; i < w; i++)
	{
		if (c[i] > cutoff)
		{
			if (ii < w2)
			for (int j = 0; j < h ; j++)
			for (int l = 0; l < pd; l++)
				y[(j*w2+ii)*pd+l] = x[(j*w+i)*pd+l];
			ii += 1;
		}
	}
}

static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i<0 || j<0 || i>=w || j>= w)
		return 0;
	else
		return x[j*w+i];
}

void reduce_width_seam_carving(float *y, int w2, float *x, int w, int h, int pd)
{
	assert(w2 > 0);
	assert(w2 < w);
	float *e = xmalloc(w * h * sizeof*e); // local energy field
	compute_energy_field(e, x, w, h, pd);
	iio_write_image_float_vec("/tmp/energy_field.npy", e, w, h, 1);

	// fill-in cumulative minimum energies
	float *m = xmalloc(w * h * sizeof*m);
	for (int i = 0; i < w*h; i++)
		m[i] = e[i];

	for (int j = 1; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float m_a = getpixel_0(m, w, h, i-1, j-1);
		float m_b = getpixel_0(m, w, h, i  , j-1);
		float m_c = getpixel_0(m, w, h, i+1, j-1);
		m[j*w+i] = e[j*w+i] + fmin(m_a, fmin(m_b, m_c));
	}
	iio_write_image_float("/tmp/carvaccum.npy", m, w, h);

}

int main(int c, char *v[])
{
	if (c < 2 || c > 4) return fprintf(stderr,
		"usage:\n\t%s columns_to_cut [in.img [out.img]]\n", *v);
		//          0 1               2       3
	int columns_to_cut = atoi(v[1]);
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	int w2 = w - columns_to_cut;
	float *y = xmalloc(w2 * h * pd * sizeof*y);

	reduce_width_naive(y, w2, x, w, h, pd);
	//reduce_width_seam_carving(y, w2, x, w, h, pd);

	iio_write_image_float_vec(filename_out, y, w2, h, pd);
	return 0;
}
