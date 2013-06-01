#include "fastmarching.h"
#include "imatge_de_ints.h"
#include "llibertat.h"

#ifndef DEBUG
#define DEBUG(...) /*fprintf(stderr,__VA_ARGS__)*/
#endif//DEBUG


//void eikonal_slow(imatge *out, imatge *in)
//{
//	error("not implemented");
//}


enum tag {KNOWN, TRIAL, FAR};

typedef struct {
	imatge *x; 		// distance function
	imatge_de_ints y[1];	// tags (known, trial, far)
	imatge_de_ints z[1];	// labels (nearest point so far)
	//int *hot;		//
	//
	int *hpos;		// locations on the heap (-42 for non-TRIAL)
	int *hvox;		// indexes of voxels on the heap
				// (filled from 0 to .nt-1)
	int nt;
} eiko_state;

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

static double dys(int a, int b, int c, int x, int y, int z)
{
	double r = hypot(x-a,hypot(y-b,z-c));
	DEBUG("\t\tdys(%d %d %d | %d %d %d) = %g\n",
			a, b, c, x, y, z, r);
	return r;
}

static double dysl(eiko_state *e, int idx, int x, int y, int z)
{
	int a = e->x->invidx[idx][0];
	int b = e->x->invidx[idx][1];
	int c = e->x->invidx[idx][2];
	return dys(a,b,c,x,y,z);
}

static double dysll(eiko_state *e, int idxa, int idxb)
{
	int a = e->x->invidx[idxa][0];
	int b = e->x->invidx[idxa][1];
	int c = e->x->invidx[idxa][2];
	return dysl(e,idxb,a,b,c);
}

static void debug_state(eiko_state *e)
{
	DEBUG("eiko state (%p):\n", (void *)e);
	DEBUG("x:");FORI(e->x->mida)DEBUG(" %g",e->x->grayplane[i]);DEBUG("\n");
	DEBUG("y:");FORI(e->y->mida)DEBUG(" %d",e->y->grayplane[i]);DEBUG("\n");
	DEBUG("z:");FORI(e->z->mida)DEBUG(" %d",e->z->grayplane[i]);DEBUG("\n");
}


#define HEAP_ENERGY(e,i) assert(i < e->nt),\
	assert(e->hvox[i] >= 0),\
	assert(e->hvox[i] < e->x->mida),\
	assert(e->y->grayplane[e->hvox[i]] == TRIAL),\
	e->x->grayplane[e->hvox[i]]

#define HEAP_SWAP(e,i,j) do{\
	int t_ = e->hvox[i];\
	e->hvox[i] = e->hvox[j];\
	e->hvox[j] = t_;\
	e->hpos[e->hvox[i]] = i;\
	e->hpos[e->hvox[j]] = j;\
}while(0)

#include "abstract_heap.h"

#define FDISPLAY_HEAP(f,e) do{\
	fprintf(f,"the heap has %d elements out of %d possible\n", e->nt, e->x->mida);\
	fprintf(f,"i\thvox\thpos\n");\
	FORI(e->nt)\
		fprintf(f,"%d\t%d\t%d\n", i, e->hvox[i], e->hpos[i]);\
	fprintf(f,"and %d off-heap elements:\n", e->x->mida - e->nt);\
	FORIR(e->nt,e->x->mida)\
		fprintf(f,"%d\t%d\t%d\n", i, e->hvox[i], e->hpos[i]);\
}while(0)
#define DISPLAY_HEAP(e) FDISPLAY_HEAP(stdout,e)

static void check_heap_consistency(eiko_state *e)
{
	//fprintf(stderr, "cheching heap consistency (with nt=%d)\n", e->nt);
	FORI(e->x->mida)
		if (e->hpos[i] != -42)
			assert(e->hvox[e->hpos[i]] == i);
	FORI(e->nt)
		assert(e->hpos[e->hvox[i]] == i);
	HEAP_ASSERT(e,e->nt);
	if(0)check_heap_consistency(e);
}

static void start_heap(eiko_state *e)
{
	XMALLOC(e->hpos, e->x->mida);
	XMALLOC(e->hvox, e->x->mida);
	FORI(e->x->mida) e->hpos[i] = -42;
	FORI(e->x->mida) e->hvox[i] = -43;
	e->nt = 0;
}

static void free_heap(eiko_state *e)
{
	xfree(e->hpos);
	xfree(e->hvox);
}

static void add_beginning_trial_to_heap(eiko_state *e, int i)
{
	//fprintf(stderr, "add_beginning_trial_to_heap %d\n", i);
	//DISPLAY_HEAP(e);
	//check_heap_consistency(e);
	assert(TRIAL == e->y->grayplane[i]);
	assert(-42 == e->hpos[i]);
	e->hvox[e->nt] = i;
	e->hpos[i] = e->nt;
	e->nt += 1;
	HEAP_ADD(e,e->nt-1);
	//check_heap_consistency(e);
	//fprintf(stderr, "\t...added!\n\n");
}

static void add_new_trial_to_heap(eiko_state *e, int i)
{
	//fprintf(stderr, "add_new_trial_to_heap %d\n", i);
	//DISPLAY_HEAP(e);
	//check_heap_consistency(e);
	assert(TRIAL == e->y->grayplane[i]);
	assert(-42 == e->hpos[i]);
	e->hvox[e->nt] = i;
	e->hpos[i] = e->nt;
	e->nt += 1;
	HEAP_ADD(e,e->nt-1);
	//check_heap_consistency(e);
	//fprintf(stderr, "\t...added!\n\n");
}


// removes it from the heap
static int pick_best_trial(eiko_state *e)
{
	//fprintf(stderr, "picking best trial...\n");
	//check_heap_consistency(e);
	assert(e->nt > 0);
	int r = e->hvox[0];
	assert(r >= 0);
	assert(r < e->x->mida);
	HEAP_REMOVE_TOP(e,e->nt);
	assert(e->hvox[e->nt-1] == r);
	e->nt -= 1;
	e->hpos[r] = -42;
	//check_heap_consistency(e);
	//fprintf(stderr, "\t...pake %d (%g)[%d]\n\n", r,
	//		e->x->grayplane[r], e->z->grayplane[r]);
	return r;
}


// ho deixa en un estat CORRECTE
static void buildup_things(eiko_state *e, imatge *x, float (*d)[3], int n)
{
	DEBUG("\n\nbuilding up things\n");
	e->x = x;
	inicialitza_himatge_de_ints(e->y, x);
	inicialitza_himatge_de_ints(e->z, x);
	x->invidx = xmalloc_invidx(x);
	e->nt = 0;
	FORI(x->mida)
	{
		e->x->grayplane[i] = INFINITY;
		e->y->grayplane[i] = FAR;
		e->z->grayplane[i] = -42;
	}
	FORI(n)
	{
		// TODO: inicialització més precisa! (sub-pixèlica)
		//float (*v)[3] = d + i;
		float *v = d[i];
		DEBUG("considering %dth point (%g %g %g)\n", i, v[0], v[1], v[2]);
		int p = v[0]; // FIXME: don't lose precision like that!
		int q = v[1]; // FIXME: don't lose precision like that!
		int r = v[2]; // FIXME: don't lose precision like that!
		if (punt_interiorP(x, p, q, r))
		{
			int idx = INDEX(x, r, q, p);
			e->x->t[r][q][p] = 0;
			e->y->t[r][q][p] = TRIAL;
			e->z->t[r][q][p] = idx;
#if 0
			int idx = INDEX(x, r, q, p);
			e->x->t[r][q][p] = 0; // FIXME: subpixel precision
			e->y->t[r][q][p] = KNOWN;
			DEBUG("IDX %d marked as KNOWN (%g)\n", idx, e->x->grayplane[idx]);
			e->z->t[r][q][p] = idx;
			int cx = 0;
			PERR(rr, r - 1, r + 2)
			PERR(qq, q - 1, q + 2)
			PERR(pp, p - 1, p + 2)
			{
				DEBUG("\t(%d %d %d) %d\n",
						pp, qq, rr, cx);
				//int idxa = INDEX(x, rr, qq, pp);
				if (punt_interiorP(x, pp, qq, rr)
					&& e->y->t[rr][qq][pp] == FAR
				   )
				{
					e->x->t[rr][qq][pp] = dysl(e,idx,
								pp,qq,rr);
					e->y->t[rr][qq][pp] = TRIAL;
					e->z->t[rr][qq][pp] = idx;
					DEBUG("\tidx %d marked trial\n", idx);
				}
				if (punt_interiorP(x, pp, qq, rr)
					&& e->y->t[rr][qq][pp] == TRIAL
				   )
				{
					// si en podem millorar la valor,
					// ho fem
					double nd = dysl(e,idx,pp,qq,rr);
					if (nd < e->x->grayplane[idx])
					{
						e->x->t[rr][qq][pp] = nd;
						e->z->t[rr][qq][pp] = idx;
					}
				}
				cx += 1;
			}
			DEBUG("cx = %d\n", cx);
			assert(cx == 27);
#endif
		}
	}
	start_heap(e);
	FORI(e->x->mida)
		if (TRIAL == e->y->grayplane[i])
			add_beginning_trial_to_heap(e, i);
}

static void free_things(eiko_state *e)
{
	allibera_imatge_de_ints(e->y);
	allibera_imatge_de_ints(e->z);
	free_heap(e);
}

#if 0
double pick_trial_distance(eiko_state *e, int idx)
{
	error("bad idea");
	assert(e->y->grayplane[idx] == TRIAL);
	double ret = INFINITY;
	int p = e->x->invidx[idx][0];
	int q = e->x->invidx[idx][1];
	int r = e->x->invidx[idx][2];
	int ivc = -42;
	PERR(rr, r - 1, r + 2)
	PERR(qq, q - 1, q + 2)
	PERR(pp, p - 1, p + 2)
	if (punt_interiorP(e->x, pp, qq, rr)
			&& e->y->t[rr][qq][pp] == KNOWN)
	{
		int i = INDEX(e->x, pp, qq, rr);
		if (ivc < 0 || e->x->grayplane[ivc] < e->x->grayplane[i])
			ivc = i;
	}
	assert(ivc >= 0);
	int pp = e->x->invidx[ivc][0];
	int qq = e->x->invidx[ivc][1];
	int rr = e->x->invidx[ivc][2];

	return ret;
}
#endif

SMART_PARAMETER(DISTANCE_SLOW,0)
void build_signed_distance(imatge *imgout, imatge *maskin)
{
	float (*points)[3] = xmalloc(maskin->mida *  sizeof(*points));
	inicialitza_himatge(imgout, maskin);
	imatge ip[2];
	FORL(2)
	{
		inicialitza_himatge(ip + l, maskin);
		//eiko_state e[1];
		int np = 0;
		FORK(maskin->d)FORJ(maskin->h)FORI(maskin->w)
			if (
				(l && maskin->t[k][j][i] > 0)
				||
				(!l && maskin->t[k][j][i] <= 0)
			)
			{
				points[np][0] = i;
				points[np][1] = j;
				points[np][2] = k;
				np += 1;
			}
		if (DISTANCE_SLOW() > 0.5)
			fill_distance_extemely_slow(ip + l, points, np);
		else
			fill_distance_fast(ip + l, points, np);
	}
	desa_imatge_general(ip, "/tmp/ip0");
	desa_imatge_general(ip+1, "/tmp/ip1");
	FORI(ip->mida)
		imgout->grayplane[i] =
		BAD_MAX(ip->grayplane[i], (ip+1)->grayplane[i]);
	desa_imatge_general(imgout, "/tmp/imgout");
	FORI(ip->mida)
		if (maskin->grayplane[i] > 0)
		       imgout->grayplane[i] *= -1;
	if (1)
	{
		FORI(ip->mida)
		{
			if (imgout->grayplane[i] > 0.6)
				imgout->grayplane[i] -= 0.5;
			if (imgout->grayplane[i] < -0.6)
				imgout->grayplane[i] += 0.5;
		}
	}
	FORL(2)
		allibera_imatge(ip + l);

	xfree(points);
}

void fill_distance_fast(imatge *x, float (*d)[3], int n)
{
	eiko_state e[1];
	buildup_things(e, x, d, n);
	DEBUG("we start with nt=%d\n", e->nt);
	debug_state(e);
	//FDISPLAY_HEAP(stderr,e);
	while (e->nt)
	{
		int i = pick_best_trial(e);
		//fprintf(stderr, "best is %d (%g) [%d]\n", i, e->x->grayplane[i], e->z->grayplane[i]);
		e->y->grayplane[i] = KNOWN;
		//fprintf(stderr, "IDX %d marked as KNOWN (%g)\n", i, e->x->grayplane[i]);
		int v[26][4], nw = fill_hinnerp_26(v, e->x, e->x->invidx[i]);
		FORK(nw)
		{
			int a = v[k][0];
			int b = v[k][1];
			int c = v[k][2];
			int q = INDEX(e->x, c, b, a);
			//if (q != v[k][3])
			//	error("%d != %d[%d %d %d]\n", q, v[k][3],
			//			v[k][0], v[k][1], v[k][2]);
			assert(q == v[k][3]);
			if (FAR == e->y->grayplane[q])
			{
				int oi = e->z->grayplane[i];
				e->x->grayplane[q] = dysl(e, oi, a, b, c);
				e->y->grayplane[q] = TRIAL;
				e->z->grayplane[q] = oi;
				add_new_trial_to_heap(e, q);
			}
			if (TRIAL == e->y->grayplane[q])
			{
				int oi = e->z->grayplane[i];
				double nv = dysll(e, oi, q);
				if (nv < e->x->grayplane[q])
				{
					//fprintf(stderr, "")
					//e->x->grayplane[q] = nv;
					e->y->grayplane[q] = TRIAL;
					e->z->grayplane[q] = oi;
					//add_trial_to_heap(e, q);
					HEAP_CHANGE_ENERGY(e,e->nt,e->hpos[q],nv);
				}
			}
		}
		//fprintf(stderr, "now nt=%d\n\n\n", e->nt);
	}
	free_things(e);
}

// x: output image
// d: array of points (must fall inside the image domain)
// n: number of points
// g: orientation vector field
void fill_hinted_signed_distance(imatge *x,
		float (*d)[3], int n, imatge *g)
{
	eiko_state e[1];
	buildup_things(e, x, d, n);
	while (e->nt)
	{
		int i = pick_best_trial(e);
		e->y->grayplane[i] = KNOWN;
		int v[26][4], nw = fill_hinnerp_26(v, e->x, e->x->invidx[i]);
		FORK(nw)
		{
			int q = v[k][3];
			if (FAR == e->y->grayplane[q])
			{
				int oi = e->z->grayplane[i];
				e->x->grayplane[q] = dysll(e, oi, q);
				e->y->grayplane[q] = TRIAL;
				e->z->grayplane[q] = oi;
				add_new_trial_to_heap(e, q);
			}
			if (TRIAL == e->y->grayplane[q])
			{
				int oi = e->z->grayplane[i];
				double nv = dysll(e, oi, q);
				if (nv < e->x->grayplane[q])
				{
					e->y->grayplane[q] = TRIAL;
					e->z->grayplane[q] = oi;
					HEAP_CHANGE_ENERGY(e,e->nt,e->hpos[q],nv);
				}
			}
		}
	}
	fprintf(stderr, "SIZES %d,  %d %d %d\n", x->mida,
			g[0].mida, g[1].mida, g[2].mida);
	FORI(3) if (x->mida != g[i].mida) error("%d != %d (%d)", x->mida, g[i].mida, i);
	FORI(x->mida)
	{
		int a = x->invidx[i][0];
		int b = x->invidx[i][1];
		int c = x->invidx[i][2];
		assert(INDEX(x,c,b,a) == i);
	}
	FORI(x->mida)
	{
		int idxp = e->z->grayplane[i];
		assert(idxp >= 0);
		assert(idxp < x->mida);
		float gp[3], vv[3];
		FORJ(3)
		{
			gp[j] = g[j].grayplane[idxp];
			vv[j] = e->x->invidx[idxp][j] - e->x->invidx[i][j];
		}
		float c = producte_escalarf(gp, vv, 3);
		float s = c < 0 ? -1 : 1;
		x->grayplane[i] *= s;
	}
	free_things(e);
}

void fill_distance_extemely_slow(imatge *x, float (*d)[3], int n)
{
	FORK(x->d)FORJ(x->h)FORI(x->w)
	{
		double od = INFINITY;
		FORL(n)
		{
			double p[3] = {i - d[l][0], j - d[l][1], k - d[l][2]};
			double nd = modul(p, 3);
			if (nd < od)
				od = nd;
		}
		x->t[k][j][i] = od;
	}
}
