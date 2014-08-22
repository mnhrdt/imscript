#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


#include "xmalloc.c"
#include "getpixel.c"

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	getsample_operator p = getsample_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < pd; l++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, pd, 2*i, 2*j, l);
		a[1] = p(in, iw, ih, pd, 2*i+1, 2*j, l);
		a[2] = p(in, iw, ih, pd, 2*i, 2*j+1, l);
		a[3] = p(in, iw, ih, pd, 2*i+1, 2*j+1, l);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[(ow*j + i)*pd+l] = cx ? m/cx : NAN;
	}
}

// evaluate a bilinear cell at the given point
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

static void bilinear_interpolation_vec(float *out, float *in,
		int w, int h, int pd,
		float p, float q)
{
	int ip = p;
	int iq = q;
	for (int l = 0; l < pd; l++) {
		float a = getsample_1(in, w, h, pd, ip  , iq  , l);
		float b = getsample_1(in, w, h, pd, ip+1, iq  , l);
		float c = getsample_1(in, w, h, pd, ip  , iq+1, l);
		float d = getsample_1(in, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p - ip, q - iq);
		out[l] = r;
	}
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		bilinear_interpolation_vec(out+pd*(ow*j+i), in, iw, ih, pd,
						(i-0.5)/2, (j-0.5)/2);
}

// y[k] = (1/n) * sum_i x[i][k]
static void float_avg(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

// euclidean distance between the vectors x and y, regularized around 0
static float fdiste(float *x, float *y, int n, float e)
{
	return n ? hypot(*x - *y, fdiste(x + 1, y + 1, n - 1, e)) : e;
}

#define WEISZ_NITER 0
#define WEISZ_EPSILON 1e-5

// y[k] = euclidean median of the vectors x[i][k]
static void float_weisz(float *y, float *x, int d, int n)
{
	float_avg(y, x, d, n);
	int niter = WEISZ_NITER;
	for (int k = 0; k < niter; k++) {
		float a[d], b = 0;
		for (int l = 0; l < d; l++)
			a[l] = 0;
		for (int i = 0; i < n; i++) {
			float dxy = fdiste(x + i*d, y, d, WEISZ_EPSILON);
			for (int l = 0; l < d; l++)
				a[l] += x[i*d + l]/dxy;
			b += 1/dxy;
		}
		for (int l = 0; l < d; l++)
			y[l] = a[l]/b;
	}
}

static void median_filter_vec(float *y, float *x, int w, int h, int pd, int rad)
{
	int np = (2 * rad + 1) * (2 * rad + 1);
	float p[np * pd];
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int cx = 0;
		for (int dy = -rad; dy <= rad; dy++)
		for (int dx = -rad; dx <= rad; dx++)
		{
			int ii = i + dx;
			int jj = j + dy;
			for (int l = 0; l < pd; l++)
				p[cx*pd+l] = getsample_1(x, w,h, pd, ii,jj, l);
			cx += 1;
		}
		float *yij = y + (j*w + i) * pd;
		float_weisz(yij, p, pd, np);
		assert(cx == np);
	}
}

static void vector_median_filter_inline(float *x, int w, int h, int pd, int rad)
{
	fprintf(stderr, "mfilter %d %d\n", w, h);
	float *tmp = xmalloc(w * h * pd * sizeof*tmp);
	median_filter_vec(tmp, x, w, h, pd, rad);
	memcpy(x, tmp, w * h * pd * sizeof*tmp);
	free(tmp);
}

static float squared_euclidean_distance(float *x, float *y, int n)
{
	int r = 0;
	for (int i = 0; i < n; i++)
	{
		float q = x[i] - y[i];
		r += q * q;
	}
	return r;
}

static float eval_displacement_by_bm(float *a, float *b, int w, int h, int pd,
		int wrad, int i, int j, float d[2])
{
	int wside = 2 * wrad + 1;
	float va[wside*wside*pd], vb[wside*wside*pd];
	int cx = 0;
	for (int dy = -wrad; dy <= wrad; dy++)
	for (int dx = -wrad; dx <= wrad; dx++)
	for (int l = 0; l < pd; l++)
	{
		int ii = i + dx;
		int jj = j + dy;
		va[cx] = getsample_0(a, w, h, pd, ii       , jj       , l);
		vb[cx] = getsample_0(b, w, h, pd, ii + d[0], jj + d[1], l);
		cx += 1;

	}
	return squared_euclidean_distance(va, vb, wside*wside*pd);
}

static void refine_displacement_at(float d[2], float *a, float *b,
		int w, int h, int pd, int wrad, int i, int j)
{
	int best_index = -1;
	int neig[13][2] = { {0,0}, {-1,-1},{-1,0},{-1,1},
		{0,-1},{0,1}, {1,-1},{1,0},{1,1}, 
		{2,0},{-2,0},{0,2},{0,-2},
	};
	float best_energy = INFINITY;

	for (int n = 0; n < 13; n++)
	{
		float D[2] = {d[0] + neig[n][0], d[1] + neig[n][1]};
		float r = eval_displacement_by_bm(a,b, w,h,pd, wrad, i,j, D);
		if (r < best_energy) {
			best_energy = r;
			best_index = n;
		}
	}
	assert(best_index >= 0);

	d[0] += neig[best_index][0];
	d[1] += neig[best_index][1];
}

static void refine_displacement(float *d, float *a, float *b,
		int w, int h, int pd, int wrad)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;
		refine_displacement_at(d + 2*idx, a, b, w, h, pd, wrad, i, j);
	}
}

void bmms_rec(float *out, float *a, float *b,
		int w, int h, int pd, int wrad, int mrad, int scale)
{
	fprintf(stderr, "scal(%d) %d %d\n", scale, w, h);
	// find an initial rhough displacement
	if (scale > 1) {
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *As = malloc(ws * hs * pd * sizeof*As);
		float *Bs = malloc(ws * hs * pd * sizeof*Bs);
		float *Os = malloc(ws * hs * 2  * sizeof*Os);
		zoom_out_by_factor_two(As, ws, hs, a, w, h, pd);
		zoom_out_by_factor_two(Bs, ws, hs, b, w, h, pd);
		bmms_rec(Os, As, Bs, ws, hs, pd, wrad, mrad, scale - 1);
		zoom_in_by_factor_two(out, w, h, Os, ws, hs, 2);
		if (mrad > 0)
			vector_median_filter_inline(out, w, h, 2, mrad);
		free(As);
		free(Bs);
		free(Os);

		for (int i = 0; i < 2*w*h; i++)
			out[i] = round(2*out[i]);
	} else {
		for (int i = 0; i < 2*w*h; i++)
			out[i] = 0;
	}

	// refine the rhough displacement by local optimization
	refine_displacement(out, a, b, w, h, pd, wrad);
}


#define MAIN_BMMS

#ifdef MAIN_BMMS
#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s WINRADIUS NSCALES MFRADIUS a.png b.png out.flo\n", *argv);
		//0 1         2       3        4     5     6
		return 1;
	}
	int winradius = atoi(argv[1]);
	int nscales = atoi(argv[2]);
	int mfradius = atoi(argv[3]);
	char *filename_a = argv[4];
	char *filename_b = argv[5];
	char *filename_out = argv[6];

	int w[2], h[2], pd[2];
	float *a = iio_read_image_float_vec(filename_a, w+0, h+0, pd+0);
	float *b = iio_read_image_float_vec(filename_b, w+1, h+1, pd+1);
	float *f = xmalloc(*w * *h * 2 * sizeof*f);

	bmms_rec(f, a, b, *w, *h, *pd, winradius, mfradius, nscales);

	iio_save_image_float_vec(filename_out, f, *w, *h, 2);

	free(a);
	free(b);
	free(f);
	return 0;
}
#endif
