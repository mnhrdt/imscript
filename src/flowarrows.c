#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "fragments.c"
#include "getpixel.c"

struct float_image {
	int w, h;
	float *x;
};

static void draw_black_pixel(int i, int j, void *data)
{
	struct float_image *x = data;
	setsample_0(x->x, x->w, x->h, 1, i, j, 0, 0);
}

static void draw_black_pixel_aa(int i, int j, float a, void *data)
{
	struct float_image *x = data;
	float v = getsample_0(x->x, x->w, x->h, 1, i, j, 0);
	v = v * (1 - a);
	setsample_0(x->x, x->w, x->h, 1, i, j, 0, v);
}

static void put_black_line(float *v, int w, int h,
		float p, float q, float r, float s)
{
	struct float_image x = {.w = w, .h = h, .x = v};
	//traverse_segment(p, q, r, s, draw_black_pixel, &x);
	traverse_segment_aa2(p, q, r, s, draw_black_pixel_aa, &x);
}

static void put_black_ball(float *v, int w, int h, float p, float q)
{
	for (int i = -1; i <= 1; i++)
	for (int j = -1; j <= 1; j++)
		setsample_0(v, w, h, 1, p+i, q+j, 0, 0);
}

#include "smapa.h"
SMART_PARAMETER(FLOWARR_MAXLEN,100)
SMART_PARAMETER(FLOWARR_MINDOT,1)
SMART_PARAMETER(FLOWARR_DODRAW,3)

static void putarrow(float *x, int w, int h, float p, float q, float u, float v)
{
	float n = hypot(u, v);
	if (n < FLOWARR_MINDOT())
		return;
	//put_black_ball(x, w, h, p, q);
	if (n > FLOWARR_MAXLEN()) {
		u *= FLOWARR_MAXLEN()/n;
		v *= FLOWARR_MAXLEN()/n;
	}
	put_black_line(x, w, h, p-u/2, q-v/2, p+u/2, q+v/2);
	if (n > FLOWARR_DODRAW()) {
		float a[2] = {p+u/2, q+v/2};
		float b[2] = {p-v/7, q+u/7};
		float c[2] = {p+v/7, q-u/7};
		put_black_line(x, w, h, a[0], a[1], b[0], b[1]);
		put_black_line(x, w, h, a[0], a[1], c[0], c[1]);
	}
}


// vv: output arrow gray image
// ff: input flow image
// s: arrow scaling
// g: grid spacing
void flowarrows(float *vv, float *ff, int w, int h, float s, int g)
{
	float (*f)[w][2] = (void*)ff;
	int gw = w/g, gh = h/g;
	for (int j = 0; j < gh; j++)
	for (int i = 0; i < gw; i++) {
		float m[2] = {0, 0}, nm = 0;
		for (int jj = 0; jj < g; jj++)
		for (int ii = 0; ii < g; ii++) {
			int pi = g*i + ii;
			int pj = g*j + jj;
			if (pi < w && pj < h && isfinite(f[pj][pi][0])
					     && isfinite(f[pj][pi][1]))	{
				m[0] += f[pj][pi][0];
				m[1] += f[pj][pi][1];
				nm += s;
			}
		}
		if (nm > 0)
			putarrow(vv, w, h, g*i+g/2, g*j+g/2, m[0]/nm, m[1]/nm);
	}
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s scale gridsize [in [out]]\n", *v);
		//                          0  1    2         3   4
		return EXIT_FAILURE;
	}
	float scale = atof(v[1]);
	int gridsize = atoi(v[2]);
	char *infile = c > 3 ? v[3] : "-";
	char *outfile = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) error("2D vector field expected");
	float *y = xmalloc(w*h*sizeof*y);
	for (int i = 0; i < w*h; i++)
		y[i] = 255;
	flowarrows(y, x, w, h, scale, gridsize);
	for (int i = 0; i   < w*h; i++)
		y[i] = (unsigned char)y[i];
	iio_write_image_float_vec(outfile, y, w, h, 1);
	return EXIT_SUCCESS;
}
#endif
