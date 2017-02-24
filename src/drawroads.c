// take a "pixelized sertit roads" data and display it as an image

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

static int insideP(int w, int h, float x, float y)
{
	return x >= 0 && y >= 0 && x < w && y < h;
}

static void render_roads(float *out, int w, int h, int offset_x, int offset_y,
		double *e, double *n, int ne, int nn)
{
	// each edge becomes a rectangle with 4 vertices
	int total_faces = ne;
	int total_vertices = 4 * ne;

	// set background to 0
	for (int i = 0; i < w*h; i++)
		out[i] = 0;

	// dump vertices
	for (int i = 0; i < ne; i++)
	{
		// a, b : integer indeces of road segment ends
		int a = e[5*i + 1];
		int b = e[5*i + 2];
		//fprintf(stderr, "ei=%d a=%d b=%d\n", i, a, b);
		assert(0 <= a); assert(a < nn);
		assert(0 <= b); assert(b < nn);
		//assert(a != b);
		assert(n[3*a + 0] == a);
		assert(n[3*b + 0] == b);

		// p, q : 2D coordinates of road segment ends
		double p[2] = { n[3*a + 1], n[3*a + 2] };
		double q[2] = { n[3*b + 1], n[3*b + 2] };
		//fprintf(stderr, "\tp = %lf %lf %lf\n", p[0], p[1], p[2]);
		//fprintf(stderr, "\tq = %lf %lf %lf\n", q[0], q[1], q[2]);

		// w : road segment width
		double rwidth = 2 * (3 - e[5*i + 4]);

		// road segment projected normal vector
		double u[2] = { q[1] - p[1], p[0] - q[0] };
		double unorm = 2 * hypot(u[0], u[1]);
		u[0] /= unorm;
		u[1] /= unorm;
		//fprintf(stderr, "\tu,w = %g %g %g\n", u[0], u[1], w);

		// v[0..3] : 3D coordinates of road segment corners
		double v[4][2] = {
			{ p[0] - rwidth * u[0], p[1] - rwidth * u[1] },
			{ p[0] + rwidth * u[0], p[1] + rwidth * u[1] },
			{ q[0] - rwidth * u[0], q[1] - rwidth * u[1] },
			{ q[0] + rwidth * u[0], q[1] + rwidth * u[1] }
		};

		// dump
		if (insideP(w, h, p[0] - offset_x, p[1] - offset_y) ||
				insideP(w, h, q[0] - offset_x, q[1] - offset_y))
		{
			printf("%g %g  %g %g  %g %g  %g %g\n",
					v[2][0] - offset_x, v[2][1] - offset_y,
					v[3][0] - offset_x, v[3][1] - offset_y,
					v[1][0] - offset_x, v[1][1] - offset_y,
					v[0][0] - offset_x, v[0][1] - offset_y
					);
		}
	}
}

#include "drawsegment.c"

struct vect_state {
	int w, h;
	float *f, vx, vy;
};

static void putvector(int i, int j, void *ee)
{
	struct vect_state *e = ee;
	if (i < 0 || j < 0 || i >= e->w || j >= e->h)
		return;
	e->f[ (i + j*e->w)*2 + 0 ] = e->vx;
	e->f[ (i + j*e->w)*2 + 1 ] = e->vy;
}

// produce a vector field of road directions
// (defined along the centre of the road, nan elsewhere)
static void render_droads(float *out, int w, int h, int offset_x, int offset_y,
		double *e, double *n, int ne, int nn)
{
	// set background to (nan,nan)
	for (int i = 0; i < 2*w*h; i++)
		out[i] = NAN;

	// dump vertices
	for (int i = 0; i < ne; i++)
	{
		// a, b : integer indeces of road segment ends
		int a = e[5*i + 1];
		int b = e[5*i + 2];
		//fprintf(stderr, "ei=%d a=%d b=%d\n", i, a, b);
		assert(0 <= a); assert(a < nn);
		assert(0 <= b); assert(b < nn);
		//assert(a != b);
		assert(n[3*a + 0] == a);
		assert(n[3*b + 0] == b);

		// p, q : 2D coordinates of road segment ends
		double p[2] = { n[3*a + 1], n[3*a + 2] };
		double q[2] = { n[3*b + 1], n[3*b + 2] };
		//fprintf(stderr, "\tp = %lf %lf %lf\n", p[0], p[1], p[2]);
		//fprintf(stderr, "\tq = %lf %lf %lf\n", q[0], q[1], q[2]);

		// v : road direction vector
		double v[2] = {q[0] - p[0], q[1] - p[1]};
		double nv = hypot(v[0], v[1]);
		v[0] /= nv;
		v[1] /= nv;

		// w : road segment width
		double rwidth = 2 * (3 - e[5*i + 4]);

		// road segment projected normal vector
		double u[2] = { q[1] - p[1], p[0] - q[0] };
		double unorm = 2 * hypot(u[0], u[1]);
		u[0] /= unorm;
		u[1] /= unorm;

		// draw segment
		int ip[2] = {round(p[0] - offset_x), round(p[1] - offset_y)};
		int iq[2] = {round(q[0] - offset_x), round(q[1] - offset_y)};
		struct vect_state vst = { w, h, out, v[0], v[1] };
		traverse_segment(ip[0], ip[1], iq[0], iq[1], putvector, &vst);
	}
}

// edges.5cols
// idx from to length2 length2
// 0   1    2  3        4
//
// nodes.3col
// idx posx posy
// 0   1    2

#include "iio.h"
#include "fail.c"
#include "parsenumbers.c"
#include "pickopt.c"
int main_drawroads(int c, char *v[])
{
	bool draw_directions = pick_option(&c, &v, "d", NULL);
	// process input arguments
	if (c != 11) {
		fprintf(stderr, "usage:\n\t%s [-d]"
		" edges.5cols nodes.3col x0 y0 w h e1 e2 e3 out.tiff\n", *v);
		//1           2          3  4  5 6 7  8  9  10
		return 1;
	}
	char *filename_e   = v[1];
	char *filename_n   = v[2];
	int param_x0       = atoi(v[3]);
	int param_y0       = atoi(v[4]);
	int w              = atoi(v[5]);
	int h              = atoi(v[6]);
	float param_e1     = atof(v[7]);
	float param_e2     = atof(v[8]);
	float param_e3     = atof(v[9]);
	char *filename_out = v[10];

	// read input lists
	int ne, nn;
	double *e = read_ascii_doubles_fn(filename_e, &ne);
	double *n = read_ascii_doubles_fn(filename_n, &nn);

	// check consistency
	if (ne%5 != 0)fail("file \"%s\" not %d cols (%d)\n", filename_e, 5, ne);
	if (nn%3 != 0)fail("file \"%s\" not %d cols (%d)\n", filename_n, 3, nn);
	ne /= 5;
	nn /= 3;

	// print debug info
	fprintf(stderr, "got %d nodes connected by %d edges\n", nn, ne);

	// render
	float *x = xmalloc(2*w*h*sizeof*x);
	if (!draw_directions) {
		render_roads(x, w, h, param_x0, param_y0, e, n, ne, nn);
		iio_write_image_float(filename_out, x, w, h);
	} else {
		render_droads(x, w, h, param_x0, param_y0, e, n, ne, nn);
		iio_write_image_float_vec(filename_out, x, w, h, 2);
	}

	// cleanup and exit
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_drawroads(c, v); }
#endif
