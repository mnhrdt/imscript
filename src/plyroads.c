// take a "sertit roads" data and display it as a ply file

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void render_roads(FILE *f, double *e, double *n, double *h,
		int ne, int nn)
{
	// each edge becomes a rectangle with 4 vertices
	int total_faces = ne;
	int total_vertices = 4 * ne;

	// ply header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by plyroads\n");
	fprintf(f, "element vertex %d\n", total_vertices);
	fprintf(f, "property double x\n");
	fprintf(f, "property double y\n");
	fprintf(f, "property double z\n");
	fprintf(f, "property uchar red\n");
	fprintf(f, "property uchar green\n");
	fprintf(f, "property uchar blue\n");
	fprintf(f, "element face %d\n", total_faces);
	fprintf(f, "property list uchar int vertex_index\n");
	fprintf(f, "end_header\n");

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

		// p, q : 3D coordinates of road segment ends
		double p[3] = { n[3*a + 1], n[3*a + 2], h[a] };
		double q[3] = { n[3*b + 1], n[3*b + 2], h[b] };
		//fprintf(stderr, "\tp = %lf %lf %lf\n", p[0], p[1], p[2]);
		//fprintf(stderr, "\tq = %lf %lf %lf\n", q[0], q[1], q[2]);

		// w : road segment width
		double w = 4 * e[5*i + 4];

		// road segment projected normal vector
		double u[2] = { q[1] - p[1], p[0] - q[0] };
		double unorm = 2 * hypot(u[0], u[1]);
		u[0] /= unorm;
		u[1] /= unorm;
		//fprintf(stderr, "\tu,w = %g %g %g\n", u[0], u[1], w);

		// v[0..3] : 3D coordinates of road segment corners
		double v[4][3] = {
			{ p[0] - w * u[0], p[1] - w * u[1], p[2] },
			{ p[0] + w * u[0], p[1] + w * u[1], p[2] },
			{ q[0] - w * u[0], q[1] - w * u[1], q[2] },
			{ q[0] + w * u[0], q[1] + w * u[1], q[2] }
		};

		// dump
		for (int j = 0; j < 4; j++)
			fprintf(f, "%lf %lf %lf 255 %d 0\n",
					v[j][0], v[j][1], v[j][2],
					(int)(100*(e[5*i + 4]-1))
					);
	}

	// output faces
	for (int i = 0; i < ne; i++)
		fprintf(f, "4 %d %d %d %d\n", 4*i, 4*i+1, 4*i+3, 4*i+2);
}

#include "fail.c"
#include "xfopen.c"
#include "parsenumbers.c"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s edges.5cols nodes.3col heights.1col o.ply", *v);
		//        0 1           2          3            4
		return 1;
	}
	char *filename_e = v[1];
	char *filename_n = v[2];
	char *filename_h = v[3];
	char *filename_out = v[4];

	// read input lists
	int ne, nn, nh;
	double *e = read_ascii_doubles_fn(filename_e, &ne);
	double *n = read_ascii_doubles_fn(filename_n, &nn);
	double *h = read_ascii_doubles_fn(filename_h, &nh);

	// check consistency
	if (ne%5 != 0)fail("file \"%s\" not %d cols (%d)\n", filename_e, 5, ne);
	if (nn%3 != 0)fail("file \"%s\" not %d cols (%d)\n", filename_n, 3, nn);
	if (nn/3 !=nh)fail("file \"%s\" not %d vals (%d)\n",filename_h,nn/3,nh);
	ne /= 5;
	nn /= 3;

	// print debug info
	fprintf(stderr, "got %d nodes connected by %d edges\n", nn, ne);

	// render
	FILE *f = xfopen(filename_out, "w");
	render_roads(f, e, n, h, ne, nn);
	xfclose(f);

	// cleanup and exit
	return 0;
}
