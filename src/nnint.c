// nearest neighbor interpolation

#include <assert.h>
#include <stdlib.h>
#include <math.h>

enum tag {KNOWN, TRIAL, FAR};

struct dist_state {
	float *x; // distance function
	int   *y; // tags (known, trial, far)
	int   *z; // labels (nearest point so far)

	int *hpos; // locations on the heap (-42 for non-TRIAL)
	int *hvox; // indexes of voxels on the heap
	           // (filled from 0 to .nt-1)
	int nt;
	int w, h;
	int (*invidx)[2];
};

#define HEAP_ENERGY(e,i) assert(i < e->nt),\
	assert(e->hvox[i] >= 0),\
	assert(e->hvox[i] < (e->w*e->h)),\
	assert(e->y[e->hvox[i]] == TRIAL),\
	e->x[e->hvox[i]]

#define HEAP_SWAP(e,i,j) do{\
	int t_ = e->hvox[i];\
	e->hvox[i] = e->hvox[j];\
	e->hvox[j] = t_;\
	e->hpos[e->hvox[i]] = i;\
	e->hpos[e->hvox[j]] = j;\
}while(0)

#include "abstract_heap.h"


// de la heap volem:
// 	* pillar el vòxel de menys valor (i saber on és)
// 	* actualitzar la valor d'un vòxel
// 	* afegir un nou vòxel amb la valor que sigui

// En un moment donat de l'algorisme,
// 	1. els vòxels "KNOWN" ja tenen els valors correctes de la distància
// 	2. els vòxels "TRIAL" tenen uns valors superiors o iguals al correcte
// 	3. els vòxels "FAR" encara no s'han recorregut
//
// en cada pas de l'algorisme:
// 	1. se selecciona el vòxel TRIAL de valor més petit
// 	2. es marca com a KNOWN i es retira de la HEAP
// 	3. s'actualitzen els vòxels TRIAL adjacents
// 	4. s'afegeien nous TRIAL adjacents

// utility function that always returns a valid pointer to memory
#include "xmalloc.c"
//static void *xmalloc(size_t n)
//{
//	void *new = malloc(n);
//	if (!new)
//	{
//		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
//		exit(1);
//	}
//	return new;
//}

static void start_heap(struct dist_state *e)
{
	e->hpos = xmalloc(e->w * e->h * sizeof*e->hpos);
	e->hvox = xmalloc(e->w * e->h * sizeof*e->hvox);
	for (int i = 0; i < e->w * e->h; i++)
	{
		e->hpos[i] = -42;
		e->hvox[i] = -43;
	}
	e->nt = 0;
}

static void free_heap(struct dist_state *e)
{
	free(e->hpos);
	free(e->hvox);
}


static void add_trial_to_heap(struct dist_state *e, int i)
{
	assert(TRIAL == e->y[i]);
	assert(-42 == e->hpos[i]);
	e->hvox[e->nt] = i;
	e->hpos[i] = e->nt;
	e->nt += 1;
	HEAP_ADD(e,e->nt-1);
}


static void fill_invidx(int (*invidx)[2], int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = w*j + i;
		invidx[idx][0] = i;
		invidx[idx][1] = j;
	}
}

// ho deixa en un estat CORRECTE
static void buildup_things(
		struct dist_state *e,
		float *dist,
		int w, int h,
		float *x
		)
{
	e->w = w;
	e->h = h;
	e->x = dist;
	e->y = xmalloc(w*h*sizeof*e->y);
	e->z = xmalloc(w*h*sizeof*e->z);
	e->invidx = xmalloc(w*h*2*sizeof(int));
	fill_invidx(e->invidx, w, h);
	e->nt = 0;
	for (int i = 0; i < w*h; i++)
	{
		e->x[i] = INFINITY;
		e->y[i] = FAR;
		e->z[i] = -42;
	}
	for (int i = 0; i < w*h; i++)
	{
		if (!isnan(x[i]))
		{
			e->x[i] = 0;
			e->y[i] = TRIAL;
			e->z[i] = i;
		}
	}
	start_heap(e);
	for (int i = 0; i < w*h; i++)
		if (TRIAL == e->y[i])
			add_trial_to_heap(e, i);
}

static void free_things(struct dist_state *e)
{
	free(e->y);
	free(e->z);
	free(e->invidx);
	free_heap(e);
}

// removes it from the heap
static int pick_best_trial(struct dist_state *e)
{
	assert(e->nt > 0);
	int r = e->hvox[0];
	assert(r >= 0);
	assert(r < e->w * e->h);
	HEAP_REMOVE_TOP(e,e->nt);
	assert(e->hvox[e->nt-1] == r);
	e->nt -= 1;
	e->hpos[r] = -42;
	return r;
}

static int get_8ineighbors(int n[8][3], int w, int h, int i, int j)
{
	int ntable[8][2] = {
		{1,0}, {0,1}, {-1,0}, {0,-1},
		{1,1}, {-1,-1}, {1,-1}, {-1,1}
	};

	int cx = 0;
	for (int p = 0; p < 8; p++)
	{
		int ii = i + ntable[p][0];
		int jj = j + ntable[p][1];
		if (ii >= 0 && jj >= 0 && ii < w && jj < h)
		{
			n[cx][0] = ii;
			n[cx][1] = jj;
			n[cx][2] = jj * w + ii;
			cx += 1;
		}
	}
	return cx;
}

static double dysl(struct dist_state *e, int idx, int x, int y)
{
	int a = e->invidx[idx][0];
	int b = e->invidx[idx][1];
	return hypot(x - a, y - b);
}

static double dysll(struct dist_state *e, int idxa, int idxb)
{
	int a = e->invidx[idxa][0];
	int b = e->invidx[idxa][1];
	return dysl(e, idxb, a, b);
}

void nnint(float *x, int w, int h)
{
	struct dist_state e[1];
	float *dist = xmalloc(w*h*sizeof*dist);
	buildup_things(e, dist, w, h, x);
	while (e->nt)
	{
		int i = pick_best_trial(e);
		e->y[i] = KNOWN;
		int v[8][3], nw = get_8ineighbors(v, w, h, e->invidx[i][0], e->invidx[i][1]);
		for (int k = 0; k < nw; k++)
		{
			int a = v[k][0];
			int b = v[k][1];
			int q = b*w + a;
			assert(q == v[k][2]);
			if (FAR == e->y[q])
			{
				int oi = e->z[i];
				e->x[q] = dysl(e, oi, a, b);
				e->y[q] = TRIAL;
				e->z[q] = oi;
				add_trial_to_heap(e, q);
			}
			if (TRIAL == e->y[q])
			{
				int oi = e->z[i];
				double nv = dysll(e, oi, q);
				if (nv < e->x[q])
				{
					e->y[q] = TRIAL;
					e->z[q] = oi;
					HEAP_CHANGE_ENERGY(e,e->nt,e->hpos[q],nv);
				}
			}
		}
	}
	for (int i = 0; i < w*h; i++)
		if (e->z[i] >= 0 && e->x[i] < w*h)
			x[i] = x[e->z[i]];
		else
			x[i] = -1;
	free_things(e);
}

void nnint_split(float *x, int w, int h, int pd)
{
	for (int l = 0; l < pd; l++)
		nnint(x + w*h*l, w, h);
}



#ifndef OMIT_NNINT_MAIN
#define USE_NNINT_MAIN
#endif

#ifdef USE_NNINT_MAIN
#include <stdio.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *filename_mask = pick_option(&c, &v, "m", "");
	int help_argument = (int)pick_option(&c, &v, "h", 0);
	if (help_argument || (c != 1 && c != 2 && c != 3)) {
		fprintf(stderr, "usage:\n\t%s [in.tiff [out.tiff]]\n", *v);
		//                          0  1        2
		return 1;
	}
	char *filename_in   = c > 1 ? v[1] : "-";
	char *filename_out  = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	if (filename_mask && *filename_mask) {
		int mw, mh;
		float *m = iio_read_image_float(filename_mask, &mw, &mh);
		for (int l = 0; l < pd; l++)
		for (int i = 0; i < mw*mh; i++)
			if (i < w*h && m[i])
				x[l*w*h+i] = NAN;
		free(m);
	}

	nnint_split(x, w, h, pd);

	iio_write_image_float_split(filename_out, x, w, h, pd);

	return 0;
}
#endif//USE_NNINT_MAIN
