// PLYROADS: take a "sertit roads" data and display it as a ply file

#include <math.h>
#include <stdio.h>

// render the SERTIT roads into a PLY file
static void render_roads(
		FILE *f,   // output file
		double *e, // array of edge data (5 columns)
		double *n, // array of node data (3 columns)
		double *h, // array of heights (1 column)
		int ne,    // number of edges
		int nn     // number of nodes
		)
{
	// each edge becomes a rectangle with 4 vertices
	int total_faces = ne;
	int total_vertices = 4 * ne;

	// print ply header
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

	// print road segments (and their colors)
	for (int i = 0; i < ne; i++)
	{
		// a, b : integer indeces of road segment ends
		int a = e[5*i + 1];
		int b = e[5*i + 2];

		// p, q : 3D coordinates of road segment ends
		double p[3] = { n[3*a + 1], n[3*a + 2], h[a] };
		double q[3] = { n[3*b + 1], n[3*b + 2], h[b] };

		// w : road segment width
		double w = 4 * e[5*i + 4];

		// road segment projected normal vector
		double u[2] = { q[1] - p[1], p[0] - q[0] };
		double unorm = 2 * hypot(u[0], u[1]);
		u[0] /= unorm;
		u[1] /= unorm;

		// v[0..3] : 3D coordinates of road segment corners
		double v[4][3] = {
			{ p[0] - w * u[0], p[1] - w * u[1], p[2] },
			{ p[0] + w * u[0], p[1] + w * u[1], p[2] },
			{ q[0] - w * u[0], q[1] - w * u[1], q[2] },
			{ q[0] + w * u[0], q[1] + w * u[1], q[2] }
		};

		// print the four coordinates
		for (int j = 0; j < 4; j++)
			fprintf(f, "%lf %lf %lf 255 %d 0\n",
				v[j][0], v[j][1], v[j][2],
				(int)(100*(e[5*i + 4]-1)) // shades of orange
				);
	}

	// print connectivity (4-faces)
	for (int i = 0; i < ne; i++)
		fprintf(f, "4 %d %d %d %d\n", 4*i, 4*i+1, 4*i+3, 4*i+2);
}

#include "parsenumbers.c"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s edges.5cols nodes.3col heights.1col o.ply\n", *v);
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

	// create output file
	FILE *f = fopen(filename_out, "w");
	if (!f) return fprintf(stderr, "ERORR: cannot write %s\n",filename_out);

	// render roads
	render_roads(f, e, n, h, ne, nn);

	// cleanup and  exit
	fclose(f);
	return 0;
}
