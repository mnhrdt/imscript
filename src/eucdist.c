// euclidean distance transform by Mejister-Felzenszwalb-Huttenlocher algorithm

#include <math.h>

static float square(float x)
{
	return x*x;
}

static float get_hypot(float *f, int *v, int q, int k)
{
	float s = (f[q] + square(q) - f[v[k]] - square(v[k])) / (2*q - 2*v[k]);
}

static void squared_distances_1d(
		float *d,  // output distance
		float *f,  // input data
		int n      // length
		)
{
	int v[n];      // index table
	float z[n+1];  // temporary array
	int k = 0;     // last index

	v[0] = 0;
	z[0] = -INFINITY;
	z[1] =  INFINITY;
	for (int q = 1; q < n; q++)
	{
		float s = get_hypot(f, v, q, k);
		while (s <= z[k])
			s = get_hypot(f, v, q, --k);
		k++;
		v[k] = q;
		z[k] = s;
		z[k+1] = INFINITY;
	}

	k = 0;
	for (int q = 0; q < n; q++)
	{
		while (z[k+1] < q)
			k++;
		d[q] = square(q - v[k]) + f[v[k]];
	}
}

static void squared_distances_2d(float *f, int w, int h)
{
	// columns
	for (int i = 0; i < w; i++)
	{
		float t[h], d[h];
		for (int j = 0; j < h; j++) t[j] = f[j*w+i];
		squared_distances_1d(d, t, h);
		for (int j = 0; j < h; j++) f[j*w+i] = d[j];
	}

	// rows
	for (int j = 0; j < h; j++)
	{
		float t[w], d[w];
		for (int i = 0; i < w; i++) t[i] = f[j*w+i];
		squared_distances_1d(d, t, h);
		for (int i = 0; i < w; i++) f[j*w+i] = d[j];
	}
}

void euclidean_distance_transform(
		float *d,   // output distance
		float *m,   // input mask (binary image)
		int w,      // width
		int h       // height
		)
{
	for (int i = 0; i < w*h; i++)
		d[i] = m[i] ? 0 : INFINITY;
	squared_distances_2d(d, w, h);
}
