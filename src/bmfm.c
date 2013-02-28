// block matching of two images guided by a fundamental matrix
// (without rectification)

void bmfm(float *disp, float *a, float *b, int w, int h, int pd, double fm[9]);


#include <math.h>

#include "xmalloc.c"
#include "getpixel.c"

static int plot_line(int (*P)[2], int w, int h, double a, double b, double c)
{
	int r = 0;
	if (fabs(a) < fabs(b)) { // slope less than 1
		double p = -a/b;
		double q = -c/b;
		for (int x = 0; x < w; x++) {
			int y = round(p*x + q);
			P[r][0] = x;
			if (y >= 0 && y < h)
				P[r++][1] = y;
		}
	} else {
		double p = -b/a;
		double q = -c/a;
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

void bmfm(float *disp, float *a, float *b, int w, int h, int pd, double fm[9])
{
	int maxpoints = 2 * (w+h);
	int (*p)[2] = xmalloc(maxpoints*sizeof*p);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (0==i&&0==j%10)fprintf(stderr,"line %d/%d\n",j,h);
		int np = plot_epipolar(p, fm, w, h, i, j);
		float mincorr = INFINITY;
		int minidx = -1;
		for (int k = 0; k < np; k++)
		{
			float c = corr(a, b, w, h, pd, i, j, p[k][0], p[k][1]);
			if (c < mincorr) {
				mincorr = c;
				minidx = k;
			}
		}
		int idx = j*w + i;
		disp[3*idx + 0] = p[minidx][0];
		disp[3*idx + 1] = p[minidx][1];
		disp[3*idx + 2] = mincorr;
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
	if (c != 5) {
		fprintf(stderr,"usage:\n\t%s a.png b.png \"fm\" out_disp\n",*v);
		//                         0 1     2       3    4
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *fm_text = v[3];
	char *filename_out = v[4];

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

	bmfm(o, a, b, w, h, pd, fm);

	iio_save_image_float_vec(filename_out, o, w, h, 3);

	free(a);
	free(b);
	free(o);

	return 0;
}

#endif//MAIN_BMFM
