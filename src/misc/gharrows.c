#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "iio.h"

#include "smapa.h"
SMART_PARAMETER(FLOWARR_MAXLEN,100)
SMART_PARAMETER(FLOWARR_MINDOT,1)
SMART_PARAMETER(FLOWARR_DODRAW,3)

static void putarrow(float p, float q, float u, float v)
{
	float n = hypot(u, v);
	if (n < FLOWARR_MINDOT())
		return;
	//put_black_ball(x, w, h, p, q);
	if (n > FLOWARR_MAXLEN()) {
		u *= FLOWARR_MAXLEN()/n;
		v *= FLOWARR_MAXLEN()/n;
	}
	printf("%g %g %g %g\n", p-u/2, q-v/2, u, v);
	//put_black_line(x, w, h, p-u/2, q-v/2, p+u/2, q+v/2);
	//if (n > FLOWARR_DODRAW()) {
	//	float a[2] = {p+u/2, q+v/2};
	//	float b[2] = {p-v/7, q+u/7};
	//	float c[2] = {p+v/7, q-u/7};
	//	put_black_line(x, w, h, a[0], a[1], b[0], b[1]);
	//	put_black_line(x, w, h, a[0], a[1], c[0], c[1]);
	//}
}



// ff: input flow image
// s: arrow scaling
// g: grid spacing
void flowarrows(float *ff, int w, int h, float s, int g)
{
	//printf("set terminal postscript eps\n");
	printf("set size ratio 1\n");
	printf("unset border\n");
	printf("unset xtics\n");
	printf("unset ytics\n");
	printf("set rmargin 0\n");
	printf("set lmargin 0\n");
	printf("set tmargin 0\n");
	printf("set bmargin 0\n");
	printf("plot [0:%d] [0:%d] \"-\" notitle with vectors "
		//	"filled "
			"linestyle 1"
			"\n", w, h);
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
			putarrow(g*i+g/2, g*j+g/2, m[0]/nm, m[1]/nm);
	}
	printf("end\n");
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	if (c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s scale gridsize [in ]\n", *v);
		//                          0  1    2         3
		return 1;
	}
	float scale = atof(v[1]);
	int gridsize = atoi(v[2]);
	char *infile = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) return fprintf(stderr,"2D vector field expected");

	flowarrows(x, w, h, scale, gridsize);



	return 0;
}
#endif
