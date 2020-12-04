// SIFT + RANSAC combined iterative matching (one iteration)

// plan:
// 	1. pick the largest 600 sift keypoints of each image
// 	2. match them exhaustively, using t=290, r=inf
// 	3. find a ransac cluster H of these pairs, n=10000, m=30, e=2
// 	4. apply H to the keypoints of the second image
// 	5. match the resulting points exhaustively, using t=290, r=2
// 	6. output the resulting list of matches



#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"

#include "siftie.c"
#define OMIT_MAIN
#include "ransac.c"
#include "ransac_cases.c"

#include "vvector.h"
static void invert_homographyf(float invH[9], float H[9])
{
	float h[3][3] = { {H[0], H[1], H[2]},
			{H[3], H[4], H[5]},
			{H[6], H[7], H[8]}};
	float det, ih[3][3];
	INVERT_3X3(ih, det, h);
	for (int i = 0; i < 9; i++) invH[i] = ih[i/3][i%3];
}

#include "smapa.h"

SMART_PARAMETER(MR_NTRIALS,10000)
SMART_PARAMETER(MR_MINLIERS,30)
SMART_PARAMETER(MR_MAXERR,2)

static int find_homographic_model_among_pairs(float *model, bool *omask,
		float *data, int n)
{
	int modeldim, datadim, nfit;
	ransac_error_evaluation_function *model_evaluation;
	ransac_model_generating_function *model_generation;
	ransac_model_accepting_function *model_acceptation = NULL;
	void *user_data = NULL;

	datadim = 4;
	modeldim = 9;
	nfit = 4;
	model_evaluation = homographic_match_error;
	model_generation = homography_from_four;
	//model_acceptation = homography_is_reasonable;
	model_acceptation = NULL;


	int ntrials = MR_NTRIALS();
	int minliers = MR_MINLIERS();
	float maxerr = MR_MAXERR();

	int n_inliers = ransac(omask, model, data, datadim, n, modeldim,
			model_evaluation, model_generation,
			nfit, ntrials, minliers, maxerr,
			model_acceptation, user_data);
	if (n_inliers > 0) {
		printf("RANSAC found a model with %d inliers\n", n_inliers);
		printf("parameters =");
		for (int i = 0; i < modeldim; i++)
			printf(" %g", model[i]);
		printf("\n");
	} else printf("GRAND FAILURE: RANSAC found no model\n");

	return n_inliers;
}

//pairs = srmatch(p[0], n[0], p[1], n[1], &npairs, t, top, rad);
static struct ann_pair *srmatch(
		int *onp,
		bool *omask,
		float oh[9],
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float t, int top, float radx, float rady)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);

	// compute exhaustive SIFT matches of the largest "top" points
	// (using a exhaustive search area)
	int tmp_na = na, tmp_nb = na, npairs;
	if (top < tmp_na) tmp_na = top;
	if (top < tmp_nb) tmp_nb = top;

	struct ann_pair *p0;
	p0 = compute_sift_matches(&npairs, ka, tmp_na, kb, tmp_nb, t,
			INFINITY, INFINITY);

	// build a temporary list with these pairs of coordinates
	float *tp0 = xmalloc(npairs*4*sizeof*tp0);
	for (int i = 0; i < npairs; i++)
	{
		struct sift_keypoint *kai = ka + p0[i].from;
		struct sift_keypoint *kbi = kb + p0[i].to;
		tp0[4*i+0] = kai->pos[0];
		tp0[4*i+1] = kai->pos[1];
		tp0[4*i+2] = kbi->pos[0];
		tp0[4*i+3] = kbi->pos[1];
	}

	// find a homographic model that matches this temporary list
	float h0[9], ih0[9];
	int nr0 = find_homographic_model_among_pairs(h0, omask, tp0, npairs);
	fprintf(stderr, "nr0 = %d\n", nr0);

	// map back all the points of the second image by this homography
	invert_homographyf(ih0, h0);
	struct sift_keypoint *hkb = xmalloc(nb*sizeof*hkb);
	memcpy(hkb, kb, nb*sizeof*hkb);
	sifthom(hkb, nb, ih0);
	//{
	//	FILE *f = xfopen("/tmp/hsifts", "w");
	//	write_raw_sifts(f, kb, nb);
	//	xfclose(f);
	//}

	// compute exhaustive SIFT matches to the registered points
	// (using a tiny search area)
	fprintf(stderr, "running the second round of matches {%g %g}...\n",
			radx, rady);
	struct ann_pair *p1;

	//p1 = compute_sift_matches(&npairs, ka, na, hkb, nb, t, radx, rady);
	int wtop = 2, htop = 2;
	for (int i = 0; i < na; i++) {
		if (ka[i].pos[0] > wtop) wtop = ka[i].pos[0];
		if (ka[i].pos[1] > htop) htop = ka[i].pos[1];
	}
	p1 = compute_sift_matches_locally(&npairs, ka, na, hkb, nb, t, radx, rady, wtop, htop);
	fprintf(stderr, "the second rund of matches found %d pairs\n", npairs);

	for (int i = 0; i < npairs; i++)
		p[i] = p1[i];

	//free(p0);

	FORI(9) oh[i] = h0[i];

	*onp = npairs;
	return p;
}

// compute pairs using sift-nn (non-sym, initial segment, explicit)
// heuristic, that starts with the first "top" keypoints to obtain an initial registration, and
int main_srmatch(int c, char *v[])
{
	if (c != 9) {
		fprintf(stderr,"usage:\n\t"
			"%s t k1 k2 top rad pairs.txt omask.txt hom.txt \n",*v);
		//        0 1 2  3  4   5   6         7         8
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[2];
	int n[2];
	FORI(2) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int npairs;
	struct ann_pair *pairs;
	float t = atof(v[1]);
	int top = atoi(v[4]);
	float rad = atof(v[5]);
	char *filename_pairs = v[6];
	char *filename_omask = v[7];
	char *filename_hom = v[8];
	float h[9] = {0};
	bool *omask = xmalloc(n[0]*sizeof*omask);
	FORI(n[0]) omask[i] = false;

	fprintf(stderr, "input read\n");
	pairs = srmatch(&npairs, omask, h, p[0],n[0], p[1],n[1], t,top,rad,rad);

	fprintf(stderr, "H = "); FORI(9) fprintf(stderr, "%g%c", h[i], (i==8)?'\n':' ');
	fprintf(stderr, "SIFTCPAIRS: produced %d pairs "
			"(from %d and %d){%d}[%g%%]\n",
			npairs, n[0], n[1], n[0]*n[1],npairs*100.0/(n[0]*n[1]));
	fprintf(stderr, "H = "); FORI(9) fprintf(stderr, "%g%c", h[i], (i==8)?'\n':' ');
	{
		FILE *f = xfopen(filename_pairs, "w");
		FORI(npairs) {
			struct sift_keypoint *ka = p[0] + pairs[i].from;
			struct sift_keypoint *kb = p[1] + pairs[i].to;
			fprintf(f, "%g %g %g %g\n",
					ka->pos[0], ka->pos[1],
					kb->pos[0], kb->pos[1]);
		}
		xfclose(f);
	}
	{
		FILE *f = xfopen(filename_omask, "w");
		FORI(npairs)
			fprintf(f, "%d%c", omask[i], (i>0&&0==i%40)?'\n':' ');
		fprintf(f, "\n");
	}
	{
		FILE *f = xfopen(filename_hom, "w");
		FORI(9) fprintf(f, "%a%c", (double)h[i], (i==8)?'\n':' ');
	}
	FORI(2) if (p[i]) free(p[i]);
	if (pairs) free(pairs);
	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_srmatch(c, v); }
#endif
