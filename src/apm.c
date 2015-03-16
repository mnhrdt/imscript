// APM, "Angulo Patch Match"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "iio.h"

// a "planet" is an image associated to an affine map
struct apm_planet {
	int w, h, pd;
	double A[3]; // affine map in the horizontal direction
	float *x;
};

// an "orbit" is a set of planets
#define MAX_ORBIT 100
struct apm_orbit {
	int n;
	struct apm_planet t[MAX_ORBIT];
};


static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	if (l < 0) l = 0;
	if (l >= pd) l = pd;
	return x[(j*w+i)*pd+l];
}

static float getsample_inf(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return INFINITY;
	if (l < 0) l = 0;
	if (l >= pd) l = pd-1;
	return x[(j*w+i)*pd+l];
}

static void horizontal_affine_warp(float *y, int yw, int yh,
		double A[3], float *x, int xw, int xh, int pd)
{
	double invA[3] = { 1/A[0], -A[1]/A[0], -A[2]/A[0] };

	// TODO: appropriate prefiltering

	for (int j = 0; j < yh; j++)
	for (int i = 0; i < yw; i++)
	for (int l = 0; l < pd; l++)
	{
		int new_i = invA[0] * i + invA[1] * j + invA[2];
		y[(j*yw+i)*pd+l] = getsample_inf(x, xw, xh, pd, new_i, j, l);
	}
}

static void find_horizontal_bounding_box(double b[2], double A[3], int w, int h)
{
	// compute the x-coordinate of each transformed corner
	double c[4] = {
		0 * A[0] + 0 * A[1] + A[2],
		w * A[0] + 0 * A[1] + A[2],
		0 * A[0] + h * A[1] + A[2],
		w * A[0] + h * A[1] + A[2]
	};

	// find the min and the max x-coordinates
	b[0] = fmin(fmin(c[0], c[1]), fmin(c[2], c[3]));
	b[1] = fmax(fmax(c[0], c[1]), fmax(c[2], c[3]));
}

static void build_planet_from_image(struct apm_planet *p,
		double A[3], float *x, int w, int h, int pd)
{
	// find required bounding box
	double bbx[2];
	find_horizontal_bounding_box(bbx, A, w, h);

	// compute image dimensions and adjust translation
	p->w = ceil(bbx[1] - bbx[0]);
	p->h = h;
	p->pd = pd;
	p->A[0] = A[0];
	p->A[1] = A[1];
	p->A[2] = A[2] - bbx[0];

	// fill-in deformed image
	p->x = malloc(p->w * p->h * pd * sizeof*p->x);
	horizontal_affine_warp(p->x, p->w, p->h, p->A, x, w, h, pd);
}

static void build_orbit_from_image(struct apm_orbit *o,
		float *x, int w, int h, int pd)
{
	double A[][3] = {
		{1, 0, 0},    // identity
		{1, -1, 0},   // negative shear, huge
		{1, 1, 0},    // positive shear, huge
		{2, 0, 0},    // compressive horizontal tilt
		{0.5, 0, 0},  // expansive horizontal tilt
		{1, 0.5, 0},  // positive shear
		{1, -0.5, 0}, // negative shear
		{1, 0.75, 0},  // positive shear
		{1, -0.75, 0}, // negative shear
		{1, 0.25, 0},  // positive shear
		{1, -0.25, 0}, // negative shear
	};
	o->n = (sizeof A)/(sizeof *A);
	assert(o->n < MAX_ORBIT);
	for (int i = 0; i < o->n; i++)
		build_planet_from_image(o->t + i, A[i], x, w, h, pd);
}

void dump_orbit_for_debugging_purposes(char *name, struct apm_orbit *o)
{
	for (int i = 0; i < o->n; i++)
	{
		char fname[FILENAME_MAX];
		snprintf(fname, FILENAME_MAX, "/tmp/apm_%s_%d.tiff", name, i);
		struct apm_planet *p = o->t + i;
		iio_save_image_float_vec(fname, p->x, p->w, p->h, p->pd);
	}
}

static double random_uniform_f(double a, double b)
{
	float u = rand()/(1.0+RAND_MAX);
	return a + (b - a) * u;
}

static int random_uniform_i(int a, int b)
{
	if (b < a)
		return random_uniform_i(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

static float eval_cost_ssd(float *a, int aw, int ah, float *b, int bw, int bh,
		int pd, int ai, int aj, int bi, int bj)
{
	double r = 0;
	int wrad = 2; // wrad=2 == 5x5 window
	for (int j = -wrad; j <= wrad; j++)
	for (int i = -wrad; i <= wrad; i++)
	for (int l = 0; l < pd; l++)
	{
		float va = getsample_inf(a, aw, ah, pd, ai+i, aj+j, l);
		if (!isfinite(va)) return INFINITY;
		float vb = getsample_inf(b, bw, bh, pd, bi+i, bj+j, l);
		if (!isfinite(vb)) return INFINITY;
		float dp = va - vb;
		r += dp * dp;
	}
	float wside = 2 * wrad + 1;
	return r/(wside*wside*wside*wside);
}

static float eval_cost(float *a, int aw, int ah, float *b, int bw, int bh,
		int pd, int ai, int aj, int bi, int bj)
{
	return eval_cost_ssd(a,aw,ah, b,bw,bh, pd, ai,aj, bi,bj);
}

static float orbital_cost(struct apm_orbit *a, struct apm_orbit *b,
		int i, int j, float disp, int aidx, int bidx)
{
	assert(aidx >= 0);
	assert(aidx < a->n);
	assert(bidx >= 0);
	assert(bidx < a->n);
	struct apm_planet *pa = a->t + aidx;
	struct apm_planet *pb = b->t + bidx;
	float *ax = pa->x; int aw = pa->w; int ah = pa->h;
	float *bx = pb->x; int bw = pb->w; int bh = pb->h;
	float ai = pa->A[0] *  i         + pa->A[1] * j + pa->A[2];
	float bi = pb->A[0] * (i + disp) + pb->A[1] * j + pb->A[2];
	int iai = lrint(ai);
	int ibi = lrint(bi);
	return eval_cost(ax,aw,ah, bx,bw,bh, pa->pd, iai,j, ibi,j);
}

#define BAD_MIN(a,b) (b)<(a)?(b):(a) 

static void init_costs_to_infinity(float *disp, float *cost,
		struct apm_orbit *a, struct apm_orbit *b)
{
	int w = a->t->w;
	int h = BAD_MIN ( a->t->h, b->t->h );

	for (int i = 0; i < w*h; i++)
		cost[i] = INFINITY;
}

static void disp_random_search(float *disp, float *cost, int *pidx,
		struct apm_orbit *a, struct apm_orbit *b,
		float *dmin, float *dmax, int ntrial)
{
	int w = a->t->w;
	int h = BAD_MIN ( a->t->h, b->t->h );

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int k = 0; k < ntrial; k++)
	{
		int ij = j * w + i;
		float d_test = random_uniform_f(dmin[ij], dmax[ij]);
		int   a_test = random_uniform_i(0, a->n - 1);
		int   b_test = random_uniform_i(0, b->n - 1);
		float c_test = orbital_cost(a, b, i, j, d_test, a_test, b_test);
		if (c_test < cost[ij]) {
			cost[ij] = c_test;
			disp[ij] = d_test;
			pidx[2*ij+0] = a_test;
			pidx[2*ij+1] = b_test;
		}
	}
}

static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

static void disp_forward_propagation(float *disp, float *cost, int *pidx,
		struct apm_orbit *a, struct apm_orbit *b,
		float *dmin, float *dmax)
{
	int w = a->t->w;
	int h = BAD_MIN ( a->t->h, b->t->h );

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j * w + i;
		//int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
		int neigs[4][2] = { {-1,0}, {-2,0}, {-1,-1}, {0,-1} };
		for (int n = 0; n < 4; n++)
		{
			int ii = i + neigs[n][0];
			int jj = j + neigs[n][1];
			if (!insideP(w, h, ii, jj)) continue;
			int nid = jj * w + ii;
			int pa = pidx[2*nid+0];
			int pb = pidx[2*nid+1];
			float d = disp[nid];
			float new_cost = orbital_cost(a, b, i, j, d, pa, pb);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				disp[idx] = d;
				pidx[2*idx+0] = pa;
				pidx[2*idx+1] = pb;
			}
		}
	}
}

static void disp_backward_propagation(float *disp, float *cost, int *pidx,
		struct apm_orbit *a, struct apm_orbit *b,
		float *dmin, float *dmax)
{
	int w = a->t->w;
	int h = BAD_MIN ( a->t->h, b->t->h ); 
	for (int i = w-1; i >= 0; i--)
	for (int j = h-1; j >= 0; j--)
	{
		int idx = j * w + i;
		//int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
		int neigs[4][2] = { {1,0}, {2,0}, {1,1}, {0,1} };
		for (int n = 0; n < 4; n++)
		{
			int ii = i + neigs[n][0];
			int jj = j + neigs[n][1];
			if (!insideP(w, h, ii, jj)) continue;
			int nid = jj * w + ii;
			int pa = pidx[2*nid+0];
			int pb = pidx[2*nid+1];
			float d = disp[nid];
			float new_cost = orbital_cost(a, b, i, j, d, pa, pb);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				disp[idx] = d;
				pidx[2*idx+0] = pa;
				pidx[2*idx+1] = pb;
			}
		}
	}
}

#include "smapa.h"
SMART_PARAMETER(APM_NITER,1)
SMART_PARAMETER(APM_NTRIAL,1)

static void orbital_apm(float *disp, float *cost, int *idx,
		struct apm_orbit *a, struct apm_orbit *b,
		float *dmin, float *dmax)
{
	int niter  = APM_NITER();
	int ntrial = APM_NTRIAL();

	init_costs_to_infinity(disp, cost, a, b);

	for (int i = 0; i < niter; i++)
	{
		disp_random_search(disp, cost, idx, a, b, dmin,dmax, ntrial);
		disp_forward_propagation(disp, cost, idx, a, b, dmin,dmax);
		disp_backward_propagation(disp, cost, idx, a, b, dmin,dmax);
	}
}

// APM : Angulo Patch Match
void apm(
		float *dout, // output disparities
		float *cout, // output costs
		int *pout,   // output orbit positions
		float *a,    // left image (reference)
		int aw,      // width
		int ah,      // height
		float *b,    // right image
		int bw,      // width
		int bh,      // height
		int pd,      // pixel dimension
		float *dmin,
		float *dmax
		// MISSING: parameters of the method:
		// ntrial, niter, winrad, comparison_method
	)
{
	// build orbits
	struct apm_orbit oa[1], ob[1];
	build_orbit_from_image(oa, a, aw, ah, pd);
	build_orbit_from_image(ob, b, bw, bh, pd);
	//dump_orbit_for_debugging_purposes("oa", oa);
	//dump_orbit_for_debugging_purposes("ob", ob);

	// run algorithm
	orbital_apm(dout, cout, pout, oa, ob, dmin, dmax);
}


int main(int c, char *v[])
{
	// process input arguments
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png dmin dmax dout.tiff cout.tiff\n", *v);
		//       0  1     2     3    4    5         6
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *str_dmin   = v[3];
	char *str_dmax   = v[4];
	char *filename_dout = v[5];
	char *filename_cout = v[6];

	// read input images
	int aw, ah, apd, bw, bh, bpd;
	float *a = iio_read_image_float_vec(filename_a, &aw, &ah, &apd);
	float *b = iio_read_image_float_vec(filename_b, &bw, &bh, &bpd);
	if (apd != bpd)
		exit(fprintf(stderr,"ERROR: image color depth mismatch\n"));

	// create dmin and dmax images
	float *dmin = malloc(aw * ah * sizeof*dmin);
	float *dmax = malloc(aw * ah * sizeof*dmax);
	for (int i = 0; i < aw * ah; i++) {
		dmin[i] = atof(str_dmin);
		dmax[i] = atof(str_dmax);
	}

	// create output images
	float *dout = malloc(aw * ah * sizeof*dout);   // disparities
	float *cout = malloc(aw * ah * sizeof*cout);   // costs
	int *pout = malloc(aw * ah * 2 * sizeof*pout); // chosen planets

	// run the algorithm
	apm(dout, cout, pout, a, aw, ah, b, bw, bh, apd, dmin, dmax);

	// save the output images
	iio_save_image_float(filename_dout, dout, aw, ah);
	iio_save_image_float(filename_cout, cout, aw, ah);
	iio_save_image_int_vec("pout.tiff", pout, aw, ah, 2);

	// cleanup and exit
	return 0;
}
