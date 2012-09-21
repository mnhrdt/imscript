#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bicubic.c"

static void cline(float *l, int n, float *x, int w, int h, float angle)
{
	float an = angle*3.1416/180;
	float c = cos(an);
	float s = sin(an);
	float t[2] = {(fmin(w,h)-1)/2.0, 0};
	float dir[2] = {c*t[0]+s*t[1], -s*t[0]+c*t[1]};
	float zer[2] = {(w-1.0)/2, (w-1.0)/2};
	float from[2], toto[2];
	for (int i = 0; i < 2; i++) {
		from[i] = zer[i] - dir[i];
		toto[i] = zer[i] + dir[i];
	}
	//fprintf(stderr, "c,s = %g %g\n", c, s);
	//fprintf(stderr, "dir = %g %g\n", dir[0], dir[1]);
	//fprintf(stderr, "zer = %g %g\n", zer[0], zer[1]);
	//fprintf(stderr, "(%g %g)=>(%g %g)\n",from[0],from[1],toto[0],toto[1]);
	for (int i = 0; i < n; i++) {
		float a = i/(n - 1.0);
		float p[2];
		for (int j = 0; j < 2; j++)
			p[j] = (1-a)*from[j] + a*toto[j];
		bicubic_interpolation(l + i, x, w, h, 1, p[0], p[1]);
	}
}

static void plot_cline(float *l, int n, char *title)
{
	if (title)
		printf("set title \"%s\"\n", title);
	printf("plot \"-\" w lines\n");
	for (int i = 0; i < n; i++)
	{
		float ipos = 2*3.1416*(i - n/2.0)/n;
		printf("\t%g %g\n", ipos, l[i]);
	}
	printf("end\n");
}

#include "iio.h"
#include "smapa.h"
SMART_PARAMETER(NFAC,1)
int main(int c, char *v[])
{
	if (c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s angle [img] >plot\n", *v);
		//                          0 1      2
		return 1;
	}
	float angle = atof(v[1]);
	char *in_img = c>2 ? v[2] : "-";
	int w, h;
	float *x = iio_read_image_float(in_img, &w, &h);
	int n = NFAC()*(2+hypot(w+2,h+2));
	float l[n];
	cline(l, n, x, w, h, angle);
	plot_cline(l, n, v[1]);
	free(x);
	return 0;
}
