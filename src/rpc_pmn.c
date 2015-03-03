// Optical flow in the altitude domain, using RPC functions
//
// The algorithm lies at the intersection of graphics, vision and remote sensing

#include <math.h>
#include <stdio.h>

#define DONT_USE_TEST_MAIN
#include "rpc.c"

#define TIFFU_OMIT_MAIN
#include "tiffu.c"

#include "xmalloc.c"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EARTH_RADIUS 6378000.0


struct ortho_view {
	// image data
	struct tiff_tile_cache *t;

	// callibration
	struct rpc *r;

	// geographic grid
	int w, h;
	double lon_0, lat_0, lon_d, lat_d;
};

// input: (i,j,h) coordinates on the raster topographic map
// output (i,j) coordinates on the image
static void project(double out[3], struct ortho_view *p, double in[3])
{
	double x = p->lon_0 + in[0] * p->lon_d;
	double y = p->lat_0 + in[1] * p->lat_d;
	assert(p->r);
	assert(out);
	eval_rpci(out, p->r, x, y, in[2]);
}

// fill-in "n" ortho views from the given point "axyh" on the first image
// The size of the geographic grid is "w,h
static void build_projection_states(
		struct ortho_view *o,      // output orthoviews
		struct tiff_tile_cache *t, // input images
		struct rpc *r,             // input rpcs
		int n,                     // number of images
		double axyh[3],            // corner on first image
		int w, int h)              // desired grid size
{
	// center in geographic coordinates
	double center[3], csouth[3];
	eval_rpc(center, r, axyh[0], axyh[1]    , axyh[2]);
	eval_rpc(csouth, r, axyh[0], axyh[1] + 1, axyh[2]);

	// stepsize given by 1 pixel to the south
	// XXX WARNING WRONG TODO FIXME : assumes North-South oriented image (!)
	double lat_step = csouth[1] - center[1];
	double latitude = center[1] * (M_PI/180);
	double lonfactor = cos(latitude);
	double lon_step = -lat_step / lonfactor;

	double lon_step_m = (lon_step*lonfactor) * (M_PI/180) * EARTH_RADIUS;
	double lat_step_m = lat_step * (M_PI/180) * EARTH_RADIUS;
	fprintf(stderr, "projection center (%g %g %g) => (%.8lf %.8lf)\n",
			axyh[0], axyh[1], axyh[2], center[0], center[1]);
	fprintf(stderr, "lon_step = %g (%g meters)\n", lon_step, lon_step_m);
	fprintf(stderr, "lat_step = %g (%g meters)\n", lat_step, lat_step_m);

	// fill-in the fields
	for (int i = 0; i < n; i++)
	{
		o[i].t = t + i;
		o[i].r = r + i;
		o[i].w = w;
		o[i].h = h;
		o[i].lon_0 = center[0];
		o[i].lat_0 = center[1];
		o[i].lon_d = lon_step;
		o[i].lat_d = lat_step;
		//o->t = t + i;
		//o->r = r + i;
		//o->w = w;
		//o->h = h;
		//o->lon_0 = center[0];
		//o->lat_0 = center[1];
		//o->lon_d = lon_step;
		//o->lat_d = lat_step;
	}
}

static double randu(double a, double b)
{
	return a + random_uniform() * (b - a);
}

#include "smapa.h"
SMART_PARAMETER(PM_WINRADIUS,1)
SMART_PARAMETER(PM_TRIALS,5)
SMART_PARAMETER(PM_NITER,3)
SMART_PARAMETER(PM_MIN,-100)
SMART_PARAMETER(PM_MAX,500)

static void huge_tiff_getpixel_float(float *out,
		struct tiff_tile_cache *t, int i, int j)
{
	void *p = tiff_tile_cache_getpixel(t, i, j);
	convert_pixel_to_float(out, t->i, p);
}

static float eval_cost_pair(
		struct tiff_tile_cache *ta, int ai, int aj,
		struct tiff_tile_cache *tb, int bi, int bj
		)
{
	int pd = ta->i->spp;

	double r = 0;
	int rad = PM_WINRADIUS();
	for (int j = -rad; j <= rad; j++)
	for (int i = -rad; i <= rad; i++)
	{
		float fa[pd]; huge_tiff_getpixel_float(fa, ta, ai + i, aj + j);
		float fb[pd]; huge_tiff_getpixel_float(fb, tb, bi + i, bj + j);
		for (int l = 0; l < pd; l++)
			r = hypot(r, fa[l] - fb[l]);
	}
	return r;
}

static float eval_cost(struct ortho_view *o, int n, int i, int j, int h)
{
	int cx = 0;
	float r = 0;
	for (int a = 0; a < n; a++)
	for (int b = 0; b < a; b++)
	{
		double ijh[3] = {i, j, h}, paijh[3], pbijh[3];
		project(paijh, o + a, ijh);
		project(pbijh, o + b, ijh);
		int ai = lrint(paijh[0]);
		int aj = lrint(paijh[1]);
		int bi = lrint(pbijh[0]);
		int bj = lrint(pbijh[1]);
		r = hypot(r, eval_cost_pair(o[a].t, ai, aj, o[b].t, bi, bj));
		cx += 1;
	}
	return r;
}
		

static void random_search(float *cost, float *out_h, float *init_h,
		struct ortho_view *t, int n)
{
	for (int k = 0; k < n; k++)
		fprintf(stderr, "rs_%d %p %p %d\n", k, t+k, t[k].r, t[k].w);

	float min_off = PM_MIN();
	float max_off = PM_MAX();
#pragma omp parallel for
	for (int j = 0; j < t->h; j++)
	for (int i = 0; i < t->w; i++)
	{
		int idx = j * t->w + i;
		for (int k = 0; k < PM_TRIALS(); k++)
		{
			float new_h = init_h[idx] + randu(min_off, max_off);
			float new_cost = eval_cost(t, n, i, j, new_h);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				out_h[idx] = new_h;
			}
		}
	}
}

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

static void backward_propagation(float *cost, float *height,
		struct ortho_view *t, int n)
{
	for (int k = 0; k < n; k++)
		fprintf(stderr, "bp_%d %p %p %d\n", k, t+k, t[k].r, t[k].w);

	int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1}};
#pragma omp parallel for
	for (int j = 0; j < t->h; j++)
	for (int i = 0; i < t->w; i++)
	{
		int idx = j * t->w + i;
		for (int k = 0; k < 4; k++)
		{
			int ii = i + neigs[k][0];
			int jj = j + neigs[k][1];
			if (!insideP(t->w, t->h, ii, jj)) continue;
			float hh = height[jj*t->w+ii];
			float new_cost = eval_cost(t, n, ii, jj, hh);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				height[idx] = hh;
			}
		}
	}
}

static void backward_propagation2(float *cost, float *height,
		struct ortho_view *t, int n)
{
	for (int k = 0; k < n; k++)
		fprintf(stderr, "bp2_%d %p %p %d\n", k, t+k, t[k].r, t[k].w);

	int neigs[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1}};
#pragma omp parallel for
	for (int i = t->w-1; i >= 0; i--)
	for (int j = t->h-1; j >= 0; j--)
	{
		int idx = j * t->w + i;
		for (int k = 0; k < 4; k++)
		{
			int ii = i + neigs[k][0];
			int jj = j + neigs[k][1];
			if (!insideP(t->w, t->h, ii, jj)) continue;
			float hh = height[jj*t->w+ii];
			float new_cost = eval_cost(t, n, ii, jj, hh);
			if (new_cost < cost[idx])
			{
				cost[idx] = new_cost;
				height[idx] = hh;
			}
		}
	}
}


#include "iio.h"
static void dump_warps(float *init_h, int w, int h,
		struct ortho_view *o, int n, char *ident)
{
	int pd = o->t->i->spp;

	for (int k = 0; k < n; k++)
	{
		struct ortho_view *ok = o + k;
		float *w0 = xmalloc(w * h * pd * sizeof(float));
		float *wh = xmalloc(w * h * pd * sizeof(float));

		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			int idx = j*w + i;
			double ij0[3] = {i, j, 0};
			double ijh[3] = {i, j, init_h[idx]};
			double p0[3], ph[3];
			project(p0, ok, ij0);
			project(ph, ok, ijh);
			int ip0[2] = {lrint(p0[0]), lrint(p0[1])};
			int iph[2] = {lrint(ph[0]), lrint(ph[1])};
			huge_tiff_getpixel_float(w0+pd*idx,ok->t,ip0[0],ip0[1]);
			huge_tiff_getpixel_float(wh+pd*idx,ok->t,iph[0],iph[1]);
		}

		char buf[FILENAME_MAX];
		snprintf(buf, FILENAME_MAX, "/tmp/pm_%s_%d_w0.tiff", ident, k);
		iio_save_image_float_vec(buf, w0, w, h, pd);
		snprintf(buf, FILENAME_MAX, "/tmp/pm_%s_%d_wh.tiff", ident, k);
		iio_save_image_float_vec(buf, wh, w, h, pd);

		free(w0);
		free(wh);
	}
}

// RPC Patch Match
void pm_rpcn(float *out_h, float *init_h, int w, int h,
		struct tiff_tile_cache *t, struct rpc *r, int n,
		double axyh[3])
{
	struct ortho_view o[n];
	build_projection_states(o, t, r, n, axyh, w, h);

	dump_warps(init_h, w, h, o, n, "X");

	float *cost = xmalloc(w * h * sizeof*cost);
	for (int i = 0; i < w * h; i++)
		cost[i] = INFINITY;

	int niter = PM_NITER();
	for (int iter = 0; iter < niter; iter++)
	{
		fprintf(stderr, "iteration %d/%d\n", iter+1, niter);
		random_search(cost, out_h, init_h, o, n);
		backward_propagation(cost, out_h , o, n);
		backward_propagation2(cost, out_h , o, n);
	}

	dump_warps(out_h, w, h, o, n, "Y");

	free(cost);
}


#define MAIN_PMN

#ifdef MAIN_PMN
#include <stdio.h>
#include "iio.h"


int main_rpc_pm(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	// input arguments
	if (c < 9 || 0 == c%2) {
		fprintf(stderr, "usage:\n\t"
		"%s (a.{tiff,rpc})+ ax ay in.tif out.tif\n", *v);
		//0  1       2     c-4 c-3 c-2    c-1
		return 1;
	}
	int n = (c - 5) / 2;
	double axyh[3] = {atof(v[c-4]), atof(v[c-3]), 0};
	char *filename_h0   = v[c-2];
	char *filename_out  = v[c-1];
	char *filename_img[n];
	char *filename_rpc[n];
	for (int i = 0; i < n; i++) {
		filename_img[i] = v[1 + 2 * i];
		filename_rpc[i] = v[2 + 2 * i];
	}

	fprintf(stderr, "patch matching from %d views\n", n);
	for (int i = 0; i < n; i++)
	{
		fprintf(stderr, "view %d/%d:\n", i+1, n);
		fprintf(stderr, "\tIMG = %s\n", filename_img[i]);
		fprintf(stderr, "\tRPC = %s\n", filename_rpc[i]);
	}

	// read input images
	int megabytes = 800/n;
	struct tiff_tile_cache t[n];
	for (int i = 0; i < n; i++)
		tiff_tile_cache_init(t + i, filename_img[i], megabytes);
	int pd = t->i->spp;
	for (int i = 0; i < n; i++)
		if (pd != t[i].i->spp)
			fail("image %d color depth mismatch\n", i);

	// read input rpcs
	struct rpc r[n];
	for (int i = 0; i < n; i++)
		read_rpc_file_xml(r + i, filename_rpc[i]);

	// read initialized raster
	int w, h;
	float *in_h0 = iio_read_image_float(filename_h0, &w, &h);

	// allocate space for output raster
	float *out_h = xmalloc(w * h * sizeof*out_h);

	// run the algorithm
	pm_rpcn(out_h, in_h0, w, h, t, r, n, axyh);

	// save the output raster
	iio_save_image_float(filename_out, out_h, w, h);

	// cleanup and exit
	free(in_h0);
	free(out_h);
	return 0;
}
int main(int c,char*v[]){return main_rpc_pm(c,v);}
#endif//MAIN_PMN
