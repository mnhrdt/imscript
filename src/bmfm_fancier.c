// block matching of two images guided by a fundamental matrix
// (without rectification)

void bmfm(float *disp, float *a, float *b, int w, int h, int pd, double fm[9]);


#include <assert.h>
#include <math.h>

#include "xmalloc.c"
#include "getpixel.c"

static int verb = 0;

static int plot_line(int (*P)[2], int w, int h, double a, double b, double c)
{
	int r = 0;
	if (fabs(a) < fabs(b)) { // slope less than 1
		double p = -a/b;
		double q = -c/b;
		if (verb) fprintf(stderr, "p q = %g %g\n", p, q);
		for (int x = 0; x < w; x++) {
			int y = round(p*x + q);
			P[r][0] = x;
			if (y >= 0 && y < h)
				P[r++][1] = y;
		}
	} else {
		double p = -b/a;
		double q = -c/a;
		if (verb) fprintf(stderr, "P Q = %g %g\n", p, q);
		for (int y = 0; y < h; y++) {
			int x = round(p*y + q);
			P[r][1] = y;
			if (x >= 0 && x < w)
				P[r++][0] = x;
		}
	}
	return r;
}

static int plot_epipolar(int (*p)[2], double fm[9], int w, int h, int i, int j)
{
	double a = i*fm[0] + j*fm[3] + fm[6];
	double b = i*fm[1] + j*fm[4] + fm[7];
	double c = i*fm[2] + j*fm[5] + fm[8];
	return plot_line(p, w, h, a, b, c);
}

static int plot_centered_segment(int (*P)[2], int w, int h,
		double a, double b, double c,
		double px, double py, double rad
		)
{
	if (!isfinite(rad)) rad = fmax(w,h);
	//fprintf(stderr, "a b c = %g %g %g\n", a, b, c);
	//fprintf(stderr, "px py = %g %g\n", px, py);
	assert(fabs(a*px + b*py + c) < 10e-6);
	double dirn = rad/hypot(a,b), sgn = a - b > 0 ? 1 : -1;
	//fprintf(stderr, "dirn sgn = %g %g\n", dirn, sgn);
	double dir[2] = {-b*dirn*sgn, a*dirn*sgn};
	//fprintf(stderr, "dir = %g %g\n", dir[0], dir[1]);
	assert(dir[0] + dir[1] > 0);
	int r = 0;
	if (fabs(a) < fabs(b)) { // slope less than 1
		double p = -a/b;
		double q = -c/b;
		int xmin = round(fmax(0,   px - dir[0]));
		int xmax = round(fmax(w-1, px + dir[0]));
		assert(xmin <= xmax);
		if (verb) fprintf(stderr, "p q = %g %g\n", p, q);
		for (int x = xmin; x <= xmax; x++) {
			int y = round(p*x + q);
			P[r][0] = x;
			if (y >= 0 && y < h)
				P[r++][1] = y;
		}
	} else {
		double p = -b/a;
		double q = -c/a;
		if (verb) fprintf(stderr, "P Q = %g %g %g %g\n", p, q, px, py);
		int ymin = round(fmax(0,   py - dir[1]));
		int ymax = round(fmin(h-1, py + dir[1]));
		//fprintf(stderr, "py dir[1] = %g %g\n", py, dir[1]);
		//fprintf(stderr, "ymin ymax = %d %d\n", ymin, ymax);
		assert(ymin <= ymax);
		for (int y = ymin; y <= ymax; y++) {
			int x = round(p*y + q);
			P[r][1] = y;
			if (x >= 0 && x < w)
				P[r++][0] = x;
		}
	}
	//fprintf(stderr, "\n");
	return r;
}

static void project_point_to_line(double pp[2],
		double a, double b, double c, double x, double y)
{
	pp[0] = (b*b*x - a*b*y - a*c)/(a*a + b*b);
	pp[1] = (a*a*y - a*b*x - b*c)/(a*a + b*b);
}

static int plot_epipolar_fancy(int (*p)[2], double fm[9], int w, int h,
		int i, int j, float ini[2], float rad)
{
	double a = i*fm[0] + j*fm[3] + fm[6];
	double b = i*fm[1] + j*fm[4] + fm[7];
	double c = i*fm[2] + j*fm[5] + fm[8];
	double ij[2] = {i, j};
	if (isfinite(*ini)) {
		ij[0] += ini[0];
		ij[1] += ini[1];
	}
	if (verb) fprintf(stderr, "a b c = %g %g %g\n", a, b, c);
	double pp[2];
	project_point_to_line(pp, a, b, c, ij[0], ij[1]);
	return plot_centered_segment(p, w, h, a, b, c, pp[0], pp[1], rad);
}

typedef float (*vector_correlation_measure)(float*,float*,int);

static float ssd_minus_mean(float *x, float *y, int n)
{
	float mx = 0, my = 0;
	for (int i = 0; i < n; i++)
	{
		mx += x[i]/n;
		my += y[i]/n;
	}

	float r = 0;
	for (int i = 0; i < n; i++)
	{
		float s = (x[i] - mx) - (y[i] - my);
		r += s * s;
	}
	return r;
}

static double corr(float *a, float *b, int w, int h, int pd,
		int ax, int ay, int bx, int by)
{
	int winside = 5;
	int n = winside*winside*pd;
	float pa[n], pb[n];
	int cx = 0;
	for (int j = 0; j < winside; j++)
	for (int i = 0; i < winside; i++)
	for (int l = 0; l < pd; l++)
	{
		int dx = i - winside/2;
		int dy = j - winside/2;
		pa[cx] = getsample_2(a, w, h, pd, ax + dx, ay + dy, l);
		pb[cx] = getsample_2(b, w, h, pd, bx + dx, by + dy, l);
		cx += 1;
	}

	vector_correlation_measure f = ssd_minus_mean;
	return f(pa, pb, n);
}

void bmfm_fancy(float *disp,         // output disparities images (dx, dy, err)
		float *a,            // input image A
		float *b,            // input image B
		int w,               // width
		int h,               // height
		int pd,              // pixel dimension
		double fm[9],        // fundamental matrix
		float *disp_init,    // initialization (optional)
		float *search_radius // optional, w.r.t. initialization
		)
{
	int maxpoints = 2 * (w+h), (*p)[2] = xmalloc(maxpoints*sizeof*p);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float rad = NAN, ini[2] = {0, 0};
		if (search_radius)
			rad = search_radius[j*w+i];
		if (disp_init) {
			ini[0] = disp_init[2*(j*w+i)+0];
			ini[1] = disp_init[2*(j*w+i)+1];
		}
		if (0==i&&0==j%10){fprintf(stderr,"line %d/%d\n",j,h);verb=1;}
		int np = plot_epipolar_fancy(p, fm, w, h, i, j, ini, rad);
		verb=0;
		float mincorr = INFINITY;
		int minidx = 0;
		for (int k = 0; k < np; k++) {
			float c = corr(a, b, w, h, pd, i, j, p[k][0], p[k][1]);
			if (c < mincorr) {
				mincorr = c;
				minidx = k;
			}
		}
		disp[3*(j*w+i) + 0] = p[minidx][0] - i;
		disp[3*(j*w+i) + 1] = p[minidx][1] - j;
		disp[3*(j*w+i) + 2] = mincorr;
	}
	free(p);
}

#ifdef MAIN_BMFM

#include <stdio.h>
#include <string.h>
#include "iio.h"

#include "fail.c"
#include "parsenumbers.c"

int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t"
				"%s a.png b.png \"fm\" out_disp rad\n", *v);
		//                0 1     2       3    4        5
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *fm_text = v[3];
	char *filename_out = v[4];
	float maxradius = atof(v[5]);

	int nfm;
	float *ffm = alloc_parse_floats(9, fm_text, &nfm);
	if (nfm != 9)
		fail("expects a fundamental matrix (9 numbers)");
	double fm[9];
	for (int i = 0; i < 9; i++)
		fm[i] = ffm[i];
	free(ffm);

	int w, h, pd, ww, hh, ppdd;
	float *a = iio_read_image_float_vec(filename_a, &w, &h, &pd);
	float *b = iio_read_image_float_vec(filename_b, &ww, &hh, &ppdd);
	if (w != ww || h != hh || pd != ppdd)
		fail("input images size mismatch");

	float *o = xmalloc(w*h*3*sizeof*o);
	float *rad = xmalloc(w*h*sizeof*o);
	for (int i = 0; i < w*h; i++)
		rad[i] = maxradius;

	bmfm_fancy(o, a, b, w, h, pd, fm, NULL, rad);

	iio_save_image_float_vec(filename_out, o, w, h, 3);

	free(a);
	free(b);
	free(o);

	return 0;
}

#endif//MAIN_BMFM
