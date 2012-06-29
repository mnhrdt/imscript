// control grid interpolation

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fail.c"
#include "xmalloc.c"
#include "bilinear_interpolation.c"
#include "random.c"



#include "smapa.h"

SMART_PARAMETER(CGI_NGRIDX,20)
SMART_PARAMETER(CGI_NGRIDY,20)
SMART_PARAMETER(CGI_MOUTITER,10)
SMART_PARAMETER(CGI_MINNITER,5)
SMART_PARAMETER(CGI_EPSILON,0.01)


struct control_point {
	float x[2]; // coordinates of the point
	float u[2]; // motion vector of that point
	int nn;     // number of neighboring points
	int n[8];   // indexes of neighboring points
	bool lopt;   // local optimality condition
};

struct control_grid {
	int ngridx, ngridy, ngrid;
	struct control_point *grid;

	// metadata: associated images
	float *f, *a, *b;
	int w, h;
};

static void print_cgi(FILE *f, struct control_grid *g)
{
	fprintf(f, "CGI %dx%d = %d points\n", g->ngridx, g->ngridy, g->ngrid);
	for (int i = 0; i < g->ngrid; i++) {
		struct control_point *p = g->grid + i;
		fprintf(f, "\tpoint %3d: (%3g %3g) {%3g %3g} [%s] %d\tns =",
				i, p->x[0], p->x[1], p->u[0], p->u[1],
				p->lopt?"true":"false", p->nn);
		for (int j = 0; j < p->nn; j++)
			fprintf(f, " %d", p->n[j]);
		fprintf(f, "\n");
	}
}

// return the top_left control point of the corresponding cell
static struct control_point *cgi_cell(struct control_grid *g, float x, float y,
	       	float *z)
{
	if (x < 0 || x > g->w-1 || y < 0 || y > g->h-1)
		return NULL;
	//fprintf(stderr, "cgi cell (%g %g)\n", x, y);
	float cellw = (g->w - 1)/(g->ngridx - 1);
	float cellh = (g->h - 1)/(g->ngridy - 1);
	int ix = x/cellw;
	int iy = y/cellh;
	int idx = iy*g->ngridx + ix;
	//fprintf(stderr, "\tix = %d\n", ix);
	//fprintf(stderr, "\tiy = %d\n", iy);
	assert(ix >= 0);
	assert(iy >= 0);
	assert(idx >= 0);
	assert(ix < g->ngridx);
	assert(iy < g->ngridy);
	assert(idx < g->ngrid);
	struct control_point *p = g->grid + idx;
	//assert(p->x[0] <= x);
	//assert(p->x[1] <= y);
	if (z) {
		z[0] = x/cellw - ix;
		z[1] = y/cellh - iy;
	}
	return p;
}

static struct control_point *cgi_neighbour(struct control_grid *g,
		struct control_point *p, int dx, int dy)
{
	int gidx = p - g->grid;
	assert(gidx >= 0);
	assert(gidx < g->ngrid);
	int x = gidx % g->ngridx;
	int y = gidx / g->ngridx;
	assert(x >= 0);
	assert(y >= 0);
	assert(x < g->ngridx);
	assert(y < g->ngridy);
	int rx = x + dx;
	int ry = y + dy;
	if (rx < 0 || rx >= g->ngridx) return NULL;
	if (ry < 0 || ry >= g->ngridy) return NULL;
	int ridx = ry * g->ngridx + rx;
	assert(ridx < g->ngrid);
	return g->grid + ridx;
}

static void cgi_eval(float *out, struct control_grid *g, float x, float y)
{
	if (x < 0 || x > g->w-1 || y < 0 || y > g->h-1) {
		out[0] = out[1] = 0;
	} else {
		float z[2];
		struct control_point *pa = cgi_cell(g, x, y, z);
		assert(pa);
		struct control_point *pb = cgi_neighbour(g, pa, 1, 0);
		struct control_point *pc = cgi_neighbour(g, pa, 0, 1);
		struct control_point *pd = cgi_neighbour(g, pa, 1, 1);
		if (!pb) pb = pa;
		if (!pc) pc = pa;
		if (!pd) pd = pa;
		for (int l = 0; l < 2; l++) {
			float a = pa->u[l];
			float b = pb->u[l];
			float c = pc->u[l];
			float d = pd->u[l];
			float r = evaluate_bilinear_cell(a, b, c, d, *z, z[1]);
			out[l] = r;
		}
	}
}

static void densify_cgi(float *f, struct control_grid *g, int w, int h)
{
	assert(w == g->w);
	assert(h == g->h);
	double eps = 0.00001;
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
			cgi_eval(f+2*(j*w+i), g, i*(1-eps), j*(1-eps));
}

static float eval_error_cell(struct control_grid *g, struct control_point *p)
{
	return 0;
}

static bool optimize_local_mv(struct control_grid *g, struct control_point *p)
{
	assert(p - g->grid >= 0);
	assert(p - g->grid < g->ngrid);
	return true;
}



void cgi(float *f, float *a, float *b, int w, int h)
{
	// build grid of control points
	struct control_grid g[1];
	g->ngridx = CGI_NGRIDX();
	g->ngridy = CGI_NGRIDY();
	g->ngrid = g->ngridx * g->ngridy;
	g->grid = xmalloc(g->ngrid * sizeof*g->grid);
	for (int j = 0; j < g->ngridy; j++)
	for (int i = 0; i < g->ngridx; i++) {
		int idx = g->ngridx * j + i;
		struct control_point *p = g->grid + idx;
		p->x[0] = i * (w - 1.0) / (g->ngridx - 1);
		p->x[1] = j * (h - 1.0) / (g->ngridy - 1);
		p->u[0] = p->u[1] = 0;
		p->nn = 0;
		for (int dj = -1; dj <= 1; dj++)
		for (int di = -1; di <= 1; di++) {
			int ii = i + di;
			int jj = j + dj;
			if (ii != i && jj != j &&
					ii >= 0 && jj >= 0 &&
					ii < g->ngridx && jj < g->ngridy)
				p->n[p->nn++] = g->ngridx * jj + ii;
		}
		p->lopt = 0;
	}
	g->f = f; g->a = a; g->b = b; g->w = w; g->h = h;

	for (int i = 0; i < g->ngrid; i++) {
		struct control_point *p = g->grid + i;
		if (p->nn == 4)
			for (int l = 0; l < 2; l++)
				p->u[l] = 1*random_normal();
	}

	print_cgi(stdout, g);

	// iterative CPMV refinement algorithm (Sullivan-Baker 1991)
	int max_iterations = CGI_MOUTITER();
	int iteration = 0;
	while(1) {
		iteration += 1;
		int count_check = 0;
		for (int j = 0; j < g->ngridy; j++)
		for (int i = 0; i < g->ngridx; i++) {
			struct control_point *p = g->grid + g->ngridx * j + i;
			if (p->lopt == 0) {
				count_check += 1;
				bool changed = optimize_local_mv(g, p);
				if (!changed)
					p->lopt = 1;
				else
					for (int k = 0; k < p->nn; k++)
						g->grid[p->n[k]].lopt = false;
			}
		}
		if (!count_check || iteration >= max_iterations) break;
	}

	densify_cgi(f, g, w, h);
	free(g->grid);
}

#ifndef OMIT_MAIN
#include "iio.h"
int main(int argc, char *argv[])
{
	if (argc != 4)
		exit(fprintf(stderr, "usage:\n\t%s a b f\n", *argv));
	char *filename_a = argv[1];
	char *filename_b = argv[2];
	char *filename_f = argv[3];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		exit(fprintf(stderr, "input images size mismatch\n"));
	float *f = xmalloc(w * h * 2 * sizeof(float));
	cgi(f, a, b, w, h);
	iio_save_image_float_vec(filename_f, f, w, h, 2);

	free(a); free(b); free(f);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
