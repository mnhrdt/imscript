#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"


#include "xmalloc.c"


static int nearest_idx(float *x, float *m, int w, int h, int pd,
		int p, int q, int width)
{
	int p0 = p - width; if (p0 < 0) p0 = 0;
	int q0 = q - width; if (q0 < 0) q0 = 0;
	int pf = p + width; if (pf >= w-1) pf = w-1;
	int qf = q + width; if (qf >= h-1) qf = h-1;
	int retval = -1;
	float dist = INFINITY;
	fprintf(stderr, "nearest idx %d %d\n", p, q);
	fprintf(stderr, "rec (%d %d) - (%d %d)\n", p0, q0, pf, qf);
	for (int j = q0; j <= qf; j++)
	for (int i = p0; i <= pf; i++)
	{
		int idx = j*w + i;
		if (m[idx] > 0 && isfinite(x[idx*pd])) {
			float newdist = hypot(i-p, j-q);
			//fprintf(stderr, "dist %d %d = %g\n", i, j, newdist);
			if (newdist < dist) {
				dist = newdist;
				retval = idx;
			}
		}
	}
	fprintf(stderr, "got %d\n", retval);
	return retval;
}

static void fill_line(float *y, float *x, int n)
{
	// compute hole locations
	int hole[n][2], nholes = 0;
	bool inhole = false;
	for (int i = 0; i < n; i++) {
		if (inhole && !isnan(x[i])) { // end of hole
			assert(i > 0 && isnan(x[i-1]));
			fprintf(stderr, "\twhile in hole %d, found a number (%g) at %d\n", nholes, x[i], i);
			hole[nholes-1][1] = i-1;
			inhole = false;
		} else if (!inhole && isnan(x[i])) { // start of hole
			if (i > 0) assert(!isnan(x[i-1]));
			fprintf(stderr, "\tout of a hole, found nan at %d, starting h=%d\n", i, nholes);
			hole[nholes][0] = i;
			inhole = true;
			nholes += 1;
		}
		y[i] = x[i];
	}
	if (inhole)
		hole[nholes-1][1] = n-1;

	fprintf(stderr, "got %d holes\n", nholes);

	// fill-in holes
	for (int h = 0; h < nholes; h++)
	{
		int a = hole[h][0];
		int b = hole[h][1];
		fprintf(stderr, "hole %d: [%d %d]\n", h, a, b);
		assert(a <= b);
		if (a > 0) assert(!isnan(x[a-1]) && isnan(x[a]));
		if (b < n-1) assert(!isnan(x[b+1]) && isnan(x[b]));
		float first = NAN, last = NAN;
		if (a > 0) first = x[a-1];
		if (b < n-1) last = x[b+1];
		if (isnan(first) && !isnan(last))
			first = last;
		if (isnan(last) && !isnan(first))
			last = first;
		fprintf(stderr, "f, l = %g, %g\n", first, last);
		if (isnan(first)) assert(isnan(last));
		if (!isnan(first)) {
			float alpha = (last - first)/(2 + b - a);
			float beta = first -  alpha * (a - 1);
			for (int i = a; i <= b; i++)
				y[i] = alpha * i + beta;
		}
	}
}

void fillcorners(float *x, float *m, int w, int h, int pd, int width)
{
	int ptab[8][2] = {
		{0,0},
		{w-1,0},
		{w-1,h-1},
		{0,h-1},
		{(w-1)/2,0},
		{(w-1)/2,h-1},
		{0,(h-1)/2},
		{w-1,(h-1)/2},
	};
	for (int pi = 0; pi < 8; pi++) {
		int p = ptab[pi][0];
		int q = ptab[pi][1];
		int nidx = nearest_idx(x, m, w, h, pd, p, q, width);
		if (nidx < 0) continue;
		fprintf(stderr, "changing (%d %d) from nidx = %d\n", p, q, nidx);
		for (int i = 0; i < pd; i++)
			x[(w*q+p)*pd + i] = x[nidx*pd + i];
		m[w*q+p] = 1;
	}

	if (pd > 1) return;

	// fill-in horizontal boundaries
	{
		float line[w];
		for (int i = 0; i < w; i++) {
			line[i] = m[i] > 0 ? x[i] : NAN;
			m[i] = 2;
		}
		fill_line(x, line, w);

		for (int i = 0; i < w; i++) {
			line[i] = m[w*(h-1)+i] > 0 ? x[w*(h-1)+i] : NAN;
			m[w*(h-1)+i] = 2;
		}
		fill_line(x+w*(h-1), line, w);
	}

	// fill-in vertical boundaries
	{
		float line[h], oline[h];
		for (int i = 0; i < h; i++) {
			line[i] = m[w*i] > 0 ? x[w*i] : NAN;
			m[w*i] = 2;
		}
		fill_line(oline, line, h);
		for (int i = 0; i < h; i++)
			x[w*i] = oline[i];

		for (int i = 0; i < h; i++) {
			line[i] = m[w*i+w-1] > 0 ? x[w*i+w-1] : NAN;
			m[w*i+w-1] = 2;
		}
		fill_line(oline, line, h);
		for (int i = 0; i < h; i++)
			x[w*i+w-1] = oline[i];
	}
}

int main(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s width i0 m0 i1 m1\n", *v);
		//                          0 1     2  3  4  5
		return EXIT_FAILURE;
	}
	int width = atoi(v[1]);
	char *in_image = v[2];
	char *in_mask = v[3];
	char *out_image = v[4];
	char *out_mask = v[5];

	int w, h, pd, ww, hh;
	float *x = iio_read_image_float_vec(in_image, &w, &h, &pd);
	float *m = iio_read_image_float(in_mask, &ww, &hh);
	if (w != ww || h != hh)
		fail("size mismatch");

	fillcorners(x, m, w, h, pd, width);

	iio_save_image_float_vec(out_image, x, w, h, pd);
	iio_save_image_float(out_mask, m, w, h);
	free(x);
	free(m);
	return EXIT_SUCCESS;
}
