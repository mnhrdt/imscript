#include <stdio.h>
#include <stdlib.h>

#include "xmalloc.c"
#include "getpixel.c"
#include "ok_list.c"
#include "grid.c"

struct ok_grid {
	struct ok_list l[1];
	struct grid g[1];
	int *buf;
};

static void ok_grid_init(struct ok_grid *o, int np,
		float x0[2], float dx[2], int n[2])
{
	int nr = n[0] * n[1];
	ok_init(o->l, nr, np);
	grid_init(o->g, 2, x0, dx, n);
	o->buf = xmalloc(np*sizeof*o->buf);
}

static void ok_grid_free(struct ok_grid *o)
{
	ok_free(o->l);
	free(o->buf);
}

static int ok_grid_add_point(struct ok_grid *o, int p, float x[2])
{
	int r = grid_locate(o->g, x);
	ok_add_point(o->l, r, p);
	return r;
}

static int ok_neighboring_points(struct ok_grid *o, float x[2])
{
	int r[4], nr = grid_locate_overlapping(r, o->g, x);
	assert(nr <= 4);
	//for (int i=0;i<nr;i++) fprintf(stderr, "\trglos{%g %g} [%d:%d] = %d\n", x[0], x[1], i, nr, r[i]);
	int cx = 0;
	for (int i = 0; i < nr; i++)
	{
		int nri = ok_which_points(o->l, r[i]);
		for (int j = 0; j < nri; j++)
			o->buf[cx++] = o->l->buf[j];
	}
	return cx;
}


static int fill_local_maxima(int *p, float *x, int w, int h, float val_min)
{
	int r = 0;
	int n[][2] = {
		{1,0}, {0,1}, {-1,0}, {0,-1}, // 4
		{1,1}, {-1,1}, {-1,-1}, {1,-1}, // 8
		{0,0}
       	}, nn = 8;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int max = 1;
		for (int k = 0; k < nn; k++)
		{
			int I = i + n[k][0];
			int J = j + n[k][1];
			float u = getpixel_1(x, w, h, i, j);
			float v = getpixel_1(x, w, h, I, J);
			if (u < val_min || u <= v)
				max = 0;
		}
		if (max)
		{
			p[2*r+0] = i;
			p[2*r+1] = j;
			r += 1;
		}
	}
	return r;
}

static void print_star_tracks(float *a, float *b, int w, int h,
		float val_min, float p_dist, float r_dist)
{
	// initialize two lists of local maxima
	int *pa = xmalloc(1*w*h*sizeof*pa);
	int *pb = xmalloc(1*w*h*sizeof*pa);
	int npa = fill_local_maxima(pa, a, w, h, val_min);
	int npb = fill_local_maxima(pb, b, w, h, val_min);
	fprintf(stderr, "found %d local maxima in image A\n", npa);
	fprintf(stderr, "found %d local maxima in image B\n", npb);

	/// build a grid structure for the second list of keypoints
	float x0[2] = {0, 0};
	float dxy[2] = {p_dist, p_dist};
	int n[2] = {1 + (w-1)/p_dist, 1 + (h-1)/p_dist};
	struct ok_grid gb[1]; ok_grid_init(gb, npb, x0, dxy, n);
	for (int i = 0; i < npb; i++)
	{
		float pos[2] = {pb[2*i+0], pb[2*i+1]};
		ok_grid_add_point(gb, i, pos);
	}

	// for each position in the first list, traverse the list of neighbors
	for (int i = 0; i < npa; i++)
	{
		float pos_a[2] = {pa[2*i+0], pa[2*i+1]};
		float val_a = a[(pa[2*i+1])*w + pa[2*i+0]];
		int nn = ok_neighboring_points(gb, pos_a);
		int *nbuf = gb->buf;

		// select closest point whose brightness is not too different
		int bestj = -1;
		float bestd = p_dist;
		for (int j = 0; j < nn; j++)
		{
			int ii = pb[2*nbuf[j]+0];
			int jj = pb[2*nbuf[j]+1];
			float pos_b[2] = {ii, jj};
			float val_b = b[jj*w + ii];
			float D = hypot(pos_b[0]-pos_a[0], pos_b[1]-pos_a[1]);
			float R = fmax(val_a/val_b, val_b/val_a);
			if (val_b > val_min && D < p_dist && R < r_dist)
			{
				bestj = j;
				bestd = D;
			}
		}
		if (bestj >= 0)
		{
			int ii = pb[2*nbuf[bestj]+0];
			int jj = pb[2*nbuf[bestj]+1];
			float pos_b[2] = {ii, jj};
			printf("%g %g %g %g\n", pos_a[0], pos_a[1],
						pos_b[0], pos_b[1]);
		}
	}

	free(pa);
	free(pb);
}

#include "iio.h"
int main(int c, char *v[])
{
	if (c != 6)
		return fprintf(stderr, "usage:\n\t"
				"%s A B minval dist valratio >pairs.txt\n", *v);
		//                0 1 2 3      4    5
	char *filename_a = v[1];
	char *filename_b = v[2];
	float minvalue = atof(v[3]);
	float distance = atof(v[4]);
	float valratio = atof(v[5]);
	int w[2], h[2];
	float *a = iio_read_image_float(filename_a, w+0, h+0);
	float *b = iio_read_image_float(filename_b, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "input size mismatch %d %d != %d %d\n",
				w[0], w[1], h[0], h[1]);
	print_star_tracks(a, b, *w, *h, minvalue, distance, valratio);
	return 0;
}
