#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fail.c"
#include "xmalloc.c"



#include "smapa.h"

SMART_PARAMETER(CGI_NGRIDX,20)
SMART_PARAMETER(CGI_NGRIDY,20)
SMART_PARAMETER(CGI_MOUTITER,10)
SMART_PARAMETER(CGI_MINNITER,5)
SMART_PARAMETER(CGI_EPSILON,0.01)


struct control_point {
	float x, y; // coordinates of the point
	float u, v; // motion vector of that point
	int nn;     // number of neighboring points
	int n[8];   // indexes of neighboring points
	bool lopt;   // local optimality condition
};

struct control_grid {
	int ngridx, ngridy, ngrid;
	struct control_point *grid;

	// associated images
	float *f, *a, *b;
	int w, h;
};

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
		p->x = i * (w - 1.0) / g->ngridx;
		p->y = j * (h - 1.0) / g->ngridy;
		p->u = p->v = 0;
		p->nn = 0;
		for (int dj = -1; dj <= 1; dj++)
		for (int di = -1; di <= 1; di++) {
			int ii = i + di;
			int jj = j + dj;
			if (ii >= 0 && jj >= 0 && ii < w && jj < h)
				p->n[p->nn++] = g->ngridx * jj + ii;
		}
		p->lopt = 0;
	}
	g->f = f; g->a = a; g->b = b; g->w = w; g->h = h;

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
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	float *f = xmalloc(w * h * 2 * sizeof(float));
	cgi(f, a, b, w, h);
	iio_save_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
