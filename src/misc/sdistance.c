void fill_distance_fast(float *dist, int w, int h, float *points, int npoints);
void fill_distance_slow(float *dist, int w, int h, float *points, int npoints);
void build_signed_distance(float *dist, float *mask, int w, int h);


#include <assert.h>
#include <stdlib.h>
#include <math.h>

void fill_distance_slow(float *dist, int w, int h, float *p, int n)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double od = INFINITY;
		for (int k = 0; k < n; k++)
		{
			double nd = hypot(i - p[2*k+0], j - p[2*k+1]);
			if (nd < od)
				od = nd;
		}
		dist[j*w+i] = od;
	}
}



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
		float *point,
		int n
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
	for (int i = 0; i < n; i++)
	{
		// TODO: inicialització més precisa! (sub-pixèlica)
		int p = point[2*i+0]; // FIXME: don't lose precision like that!
		int q = point[2*i+1]; // FIXME: don't lose precision like that!
		if (p >= 0 && p < w && q >= 0 && q < h)
		{
			int idx = w*q + p;
			e->x[idx] = 0;
			e->y[idx] = TRIAL;
			e->z[idx] = idx;
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

void fill_distance_fast(
		float *dist, // output image, to be filled with distances
		int w,       // width of ouput image
		int h,       // height of input image
		float *p,    // list of input point coordinates
		int n        // number of  input points
		)
{
	struct dist_state e[1];
	buildup_things(e, dist, w, h, p, n);
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
	free_things(e);
}


#include "iio.h"
void build_signed_distance(float *d, float *m, int w, int h)
{
	iio_write_image_float("/tmp/mmmm.tiff", m, w, h);
	float *t = xmalloc(    w * h * sizeof*t);  // temporary image
	float *p = xmalloc(2 * w * h * sizeof*p);  // list of points
	for (int l = 0; l < 2; l++)
	{
		int n = 0;
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			if ( l && !m[j*w+i]) continue;
			if (!l &&  m[j*w+i]) continue;
			p[2*n+0] = i;
			p[2*n+1] = j;
			n += 1;
		}
		fill_distance_fast(t, w, h, p, n);
		if ( l) iio_write_image_float("/tmp/sdistemp_1.tiff", t, w, h);
		if (!l) iio_write_image_float("/tmp/sdistemp_0.tiff", t, w, h);
		for (int i = 0; i < w*h; i++)
		{
			if ( l && !m[i]) d[i] = t[i];
			if (!l &&  m[i]) d[i] = -t[i];
		}
	}
	free(p);
	free(t);
}



#ifndef OMIT_DISTANCE_MAIN
#define USE_DISTANCE_MAIN
#endif


#ifdef USE_DISTANCE_MAIN
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main_sdistance(int c, char *v[])
{
	if (c > 1)
		return fprintf(stderr,
				"usage:\n\t%s <masrk >signed_distance\n", *v);
		//                          0

	int w, h;
	float *m = iio_read_image_float("-", &w, &h);
	float *d = malloc(w*h*sizeof*d);
	build_signed_distance(d, m, w, h);
	iio_write_image_float("-", d, w, h);
	return 0;
}
#endif//USE_DISTANCE_MAIN

#ifndef OMIT_ALL_MAINS
int main(int c, char *v[]) { return main_sdistance(c-1, v+1); }
#endif//OMIT_ALL_MAINS
