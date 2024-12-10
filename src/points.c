// various operations with point clouds (generation, filtering)

// TODO: complete with more n-dimensional configurations
// * uniform points on the surface of a sphere
// * uniform points on the interior of a sphere (hard!)
// * apply a random element of SO(d)
// * apply a random element of SL(d)
// * apply a random element of GL(d)
// * U*Σ with prescribed diagnonal of Σ
// * N points uniformly _spaced_ on the edges of
//   - the d-simplex
//   - the d-hypercube
//   - the d-cross-polytope
//   - the coordinate axes
//
// * Possible in-band notation for "edges": a NAN point means that the next 2
// points are an edge.  A point of the form (NAN, k, ...), means that the next
// k points form a k-simplex.

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"
#include "drawsegment.c"
#include "pickopt.c"

#define π 3.14159265358979323846264338328


#include "random.c"



// struct tag to hold the data for a multidimensional wireframe
struct wireframe {
	int d;     // dimension
	int n;     // number of points
	int m;     // number of edges
	float *x;  // coordinates of the points (d*n floats)
	int *e;  // indices of the edges (2*m ints)
};

static void fprint_wireframe(FILE *f, struct wireframe *w)
{
	fprintf(f, "wireframe\n");
	fprintf(f, "d = %d\n", w->d);
	fprintf(f, "n = %d\n", w->n);
	fprintf(f, "m = %d\n", w->m);
	fprintf(f, "x =\n");
	for (int j = 0; j < w->n; j++)
	{
		for (int i = 0; i < w->d; i++)
			fprintf(f, "\t%g", w->x[j*w->d+i]);
		fprintf(f, "\n");
	}
	for (int j = 0; j < w->m; j++)
		fprintf(f, "\t%d\t%d\n", w->e[2*j+0], w->e[2*j+1]);
}

// replace each edge of a wireframe with s>=3 uniformly distributed points
static void unwireframize(struct wireframe *w, int s)
{
	assert(s >= 3);

	// the original points are kept, for each edge we add s-2 points inside
	int n = w->n + (s - 2) * w->m;
	float *x = xmalloc(n * w->d * sizeof*x);
	memcpy(x, w->x, w->n * w->d * sizeof*x);
	for (int j = 0; j < w->m; j++)
	{
		int a = w->e[2*j+0];
		int b = w->e[2*j+1];
		float *p = x + a * w->d;
		float *q = x + b * w->d;
		for (int k = 1; k < s - 1; k++)
		{
			float t = k/(s-1.0);
			float *y = x + w->n * w->d;
			for (int i = 0; i < w->d; i++)
				y[i] = t*p[i] + (1 - t)*q[i];
			w->n += 1;
		}
	}

	// remove previous edges
	xfree(w->e);
	w->m = 0;
	w->e = NULL;

	// record new points
	xfree(w->x);
	w->n = n;
	w->x = x;
}

static void build_positive_axes(struct wireframe *w, int d)
{
	w->d = d;         // dimension
	w->n = d + 1;     // one point per dimension, plus the center 0
	w->m = d;         // one edge per dimension
	w->x = xmalloc(w->n * d * sizeof*w->x);
	w->e = xmalloc(w->m * 2 * sizeof*w->e);
	for (int i = 0; i < w->n * d; i++) // zero everything
		w->x[i] = 0;
	for (int j = 0; j < d; j++)        // jth component of x_j is 1
		w->x[d*j + j] = 1;
	for (int j = 0; j < d; j++)        // edges from 0 to each e_j
	{
		w->e[2*j+0] = d;           // index of 0
		w->e[2*j+1] = j;           // index of e_j
	}
}

static void build_centered_axes(struct wireframe *w, int d)
{
	w->d = d;          // dimension
	w->n = 2*d + 1;    // two points per axis, plus the 0
	w->m = 2*d;        // two edges per axis
	w->x = xmalloc(w->n * d * sizeof*w->x);
	w->e = xmalloc(w->m * 2 * sizeof*w->e);
	for (int i = 0; i < w->n * d; i++)  // zero all
		w->x[i] = 0;
	for (int j = 0; j < d; j++)
	{
		w->x[d*(2*j+0) + j] = 1;
		w->x[d*(2*j+1) + j] = -1;
	}
	for (int j = 0; j < d; j++)
	{
		w->e[2*(2*j+0)+0] = 2*d;
		w->e[2*(2*j+0)+1] = 2*j + 0;
		w->e[2*(2*j+1)+0] = 2*d;
		w->e[2*(2*j+1)+1] = 2*j + 1;
	}
}

// TODO: add builders for all regular polytopes in dimensions 4 and above

static void print_points(float *x, int n, int d)
{
	for (int i = 0; i < n; i++)
	for (int j = 0; j < d; j++)
		printf("%g%c", x[i*d+j], (j==d-1)?'\n':' ');
}

#include "smapa.h"
SMART_PARAMETER_SILENT(SRAND,0)

int main_random(int c, char *v[])
{
	float s = atof(pick_option(&c, &v, "s", "1"));
	char *offset_string = pick_option(&c, &v, "o", "");
	//fprintf(stderr, "c=%d\n", c);
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s dimension distribution npoints"
		//                          0 1         2            3
		" [-s param] [-t translate]\n", *v);
		return 1;
	}
	int d = atoi(v[1]);
	int n = atoi(v[3]);
	//fprintf(stderr, "generating %d %d-dimensional points from \"%s\"\n",
	//		n, d, v[2]);
	float *x = xmalloc(d*n*sizeof*x);
	double offset[d];
	for (int i = 0; i < d; i++)
		offset[i] = 0;
	if (*offset_string)
		read_n_doubles_from_string(offset, offset_string, d);

	xsrand(100+SRAND());

	switch (v[2][0]) {
	case 'g':
		for (int i = 0; i < n*d; i++)
			x[i] = s * random_normal();
		break;
	case 'u':
		for (int i = 0; i < n*d; i++)
			x[i] = s * (random_uniform()-0.5);
		break;
	case 'c':
		if (d != 2) fail("%d-dimensional cauchy not implemented\n", d);
		for (int i = 0; i < n; i++)
		{
			double θ = 2 * π * random_uniform();
			double ρ = s * random_cauchy();
			x[2*i+0] = ρ * cos(θ);
			x[2*i+1] = ρ * sin(θ);
		}
		break;
	default:
		return fprintf(stderr,"unrecognized distribution \"%s\"",v[2]);
	}

	for (int i = 0; i < n; i++)
	for (int j = 0; j < d; j++)
		x[i*d+j] += offset[j];

	print_points(x, n, d);

	return 0;
}

// map a set of points (by default, the identity)
int main_map(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr,"usage:\n\t%s type params <in >out\n",*v);
		//                                0 1    2
	char *map_type   = v[1];
	char *map_params = v[2];

	int nparams;
	double *params = alloc_parse_doubles(10000, map_params, &nparams);

	int n, d;
	float *x = iio_read_image_float("-", &d, &n);

	switch (*map_type) {
	case 't': // translation
		if (nparams != d) fail("translation d=%d np=%d", d, nparams);
		for (int i = 0; i < n; i++)
		for (int j = 0; j < d; j++)
			x[i*d+j] += params[j];
		break;
	case 's': // scaling
		if (nparams != 1) fail("scaling d=%d np=%d", d, nparams);
		for (int i = 0; i < n*d; i++)
			x[i] *= params[0];
		break;
	default:
		fail("unrecognized map type \"%s\"", map_type);
	}


	print_points(x, n, d);

	return 0;
}

int main_polytope(int c, char *v[])
{
	// named arguments
	int s = atoi(pick_option(&c, &v, "s", "10"));
	if (!s) return fprintf(stderr, "please put an s, not %d\n", s);

	// positional arguments
	if (c != 3)
		return fprintf(stderr, "usage\n\t%s dim type\n", *v);
		//                                0 1   2

	int d = atoi(v[1]); // dimension
	char *poly_type = v[2];

	struct wireframe p[1];
	//build_positive_axes(p, d);
	build_centered_axes(p, d);
	//fprint_wireframe(stderr, p);
	unwireframize(p, s);
	print_points(p->x, p->n, d);

	return 0;
}

// special configurations
int main_config(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr,"usage:\n\t%s type params\n",*v);
		//                                0 1    2
	char *config_type   = v[1];
	char *config_params = v[2];

	int nparams;
	double *p = alloc_parse_doubles(10000, config_params, &nparams);

	int n = -1;
	int d = -1;

	if (false) { ;
	} else if (0 == strcmp(config_type, "twosquares")) {
		// two squares in the "difficult" configuration
		if (nparams != 2) fail("bad twosquares np=%d", nparams);
		if (isnan(p[0])) p[0] = 10;
		if (isnan(p[1])) p[1] = 0.25;
		n = 8;
		d = 2;
		float x[8][2] = {
			{0,0}, {0,1}, {1,0}, {1,1},
			{p[0]       , p[0]       },
			{p[0] + p[1], p[0]       },
			{p[0]       , p[0] + p[1]},
			{p[0] + p[1], p[0] + p[1]},

		};
		print_points(*x, n, d);
	} else if (0 == strcmp(config_type, "triangle")) {
		// a triangle of side 1, angle p[0] in degrees, side p[2]
		if (nparams != 2) fail("bad triangle np=%d", nparams);
		if (isnan(p[0])) p[0] = 90;
		if (isnan(p[1])) p[1] = 1;
		n = 3;
		d = 2;
		float x[3][2] = {
			{0, 0}, {1, 0},
			{p[1]*cos(π*p[0]/180), p[1]*sin(π*p[0]/180)},
		};
		print_points(*x, n, d);
	} else
		fail("unrecognized config \"\"", config_type);

	return 0;
}



// CLI utility to access some point processing programs
int main_points(int c, char *v[])
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "random")) return main_random(c-1, v+1);
	else if (0 == strcmp(v[1], "map")) return main_map(c-1, v+1);
	else if (0 == strcmp(v[1], "config")) return main_config(c-1, v+1);
	else if (0 == strcmp(v[1], "polytope")) return main_polytope(c-1, v+1);
//	else if (0 == strcmp(v[1], "stats")) return main_stats(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s [random|config|stats] "
			       "params... \n", *v);
	       return 1;
	}
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_points(c, v); }
#endif
