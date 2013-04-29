// palette: a program to tranform from GRAY FLOATS to RGB BYTES
// SYNOPSIS:
// 	palette FROM TO PALSPEC [in.gray [out.rgb]]
//
//	FROM    = real value (-inf=min)
//	TO      = real value (inf=max)
//	PALSPEC = named palette (HOT, GLOBE, RELIEF, etc)j

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>


//#define PALSAMPLES 1024
#define PALSAMPLES 256

// Idea: The function
//
// 	palette : R -> rgb
//
// is described as a composition of two steps
//
// 	quantization : R -> {0, ..., PALSAMPLES-1}
// 	lookup : {0, ..., PALSAMPLES} -> rgb
//
// The first step is an affine function followed by a rounding, and the second
// step is the evaluation of a look-up table.  The proposed implementation is
// kept simple and the size of the look-up table is fixed.




#include "fail.c"

struct palette {
	uint8_t t[3*PALSAMPLES];
	float m, M;
};

static void fprint_palette(FILE *f, struct palette *p)
{
	fprintf(f, "palette %g %g\n", p->m, p->M);
	for (int i = 0; i < PALSAMPLES; i++)
		fprintf(f, "\tpal[%d] = {%d %d %d}\n",
			i, p->t[3*i+0], p->t[3*i+1], p->t[3*i+2]);
}

static void fill_palette_gray(uint8_t *t)
{
	float factor = 255.0 / PALSAMPLES;
	for (int i = 0; i < PALSAMPLES; i++)
	for (int j = 0; j < 3; j++)
		t[3*i+j] = factor * i;
}

static void fill_palette_hot(uint8_t *t)
{
	float factor = 255.0 / PALSAMPLES;
	for (int i = 0; i < PALSAMPLES; i++)
	{
		t[3*i+0] = 255;
		if (i < 85)
			t[3*i+0] = 3 * factor * i;
		t[3*i+1] = factor * i;
		t[3*i+2] = factor * i;
	}
}

static void fill_palette_with_nodes(struct palette *p, float *n, int nn)
{
	if (nn < 2) fail("at least 2 nodes required");
	float factor = (PALSAMPLES-1.0)/(nn-1.0);
	assert(PALSAMPLES == 1+factor*(nn-1));
	for (int i = 0; i < nn-1; i++)
	{
		float posi = n[4*i];
		float posnext = n[4*(i+1)];
		fprintf(stderr, "fwpn node %d : %g %g\n", i, posi, posnext);
		if (posi > posnext)
			fail("palette nodes must be in increasing order");
		int idxi = factor * i;
		int idxnext = factor * (i+1);
		assert(idxi <= idxnext);
		assert(idxnext < PALSAMPLES);
		for (int j = 0; j < (idxnext-idxi); j++)
		{
			float a = j*1.0/(idxnext - idxi);
			p->t[3*(idxi+j)+0] = (1-a)*n[4*i+1] + a*n[4*(i+1)+1];
			p->t[3*(idxi+j)+1] = (1-a)*n[4*i+2] + a*n[4*(i+1)+2];
			p->t[3*(idxi+j)+2] = (1-a)*n[4*i+3] + a*n[4*(i+1)+3];
		}
	}
	p->m = n[0];
	p->M = n[4*(nn-1)];
}

static void set_node_positions_linearly(float *n, int nn, float m, float M)
{
	float beta = (M - m)/(nn - 1);
	for (int i = 0; i < nn; i++)
		n[4*i] = m + beta * i;
}

static float nodes_cocoterrain[] = {
	1,   0,   0,   0, // black
	2, 255,   0, 255, // magenta
	3,   0,   0, 255, // blue
	4,   0, 255, 255, // cyan
	5,   0, 255,   0, // green
	6, 170,  84,   0, // brown
	7, 255, 255, 255, // white
};


static void fill_palette(struct palette *p, char *s, float m, float M)
{
	p->m = m;
	p->M = M;
	if (0 == strcmp(s, "gray"))
		fill_palette_gray(p->t);
	else if (0 == strcmp(s, "hot"))
		fill_palette_hot(p->t);
	else if (0 == strcmp(s, "cocoterrain")) {
		set_node_positions_linearly(nodes_cocoterrain, 7, m, M);
		fill_palette_with_nodes(p, nodes_cocoterrain, 7);
	}
	else fail("unrecognized palette \"%s\"", s);
}

static void get_palette_color(uint8_t *rgb, struct palette *p, float x)
{
	int ix = round((PALSAMPLES-1)*(x - p->m)/(p->M - p->m));
	if (ix < 0) ix = 0;
	if (ix >= PALSAMPLES) ix = PALSAMPLES - 1;
	rgb[0] = p->t[3*ix+0];
	rgb[1] = p->t[3*ix+1];
	rgb[2] = p->t[3*ix+2];
}

static void get_min_max(float *min, float *max, float *x, int n)
{
	float m = INFINITY, M = -m;
	for (int i = 0; i < n; i++)
		if (isfinite(x[i]))
		{
			m = fmin(m, x[i]);
			M = fmax(M, x[i]);
		}
	if (min) *min = m;
	if (max) *max = M;
}

void apply_palette(uint8_t *y, float *x, int n, char *s, float m, float M)
{
	if (!isfinite(m)) get_min_max(&m, 0, x, n);
	if (!isfinite(M)) get_min_max(0, &M, x, n);

	struct palette p[1];
	fill_palette(p, s, m, M);

	fprint_palette(stderr, p);

	for (int i = 0; i < n; i++)
		get_palette_color(y + 3*i, p, x[i]);
}

#define PALETTE_MAIN

#ifdef PALETTE_MAIN
#include "iio.h"
#include "xmalloc.c"

int main(int c, char *v[])
{
	if (c != 4 && c != 5 && c != 6 ) {
		fprintf(stderr, "usage:\n\t%s from to pal [in [out]]\n", *v);
		//                         0  1    2   3    4   5
		return 1;
	}
	float from = atof(v[1]);
	float to = atof(v[2]);
	char *palette_id = v[3];
	char *filename_in = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	uint8_t *y = xmalloc(3*w*h);

	apply_palette(y, x, w*h, palette_id, from, to);

	iio_save_image_uint8_vec(filename_out, y, w, h, 3);

	free(x);
	free(y);
	return 0;
}
#endif//PALETTE_MAIN
