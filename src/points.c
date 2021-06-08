// various operations with point clouds (generation, filtering)

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"
#include "drawsegment.c"
#include "pickopt.c"

#define π 3.14159265358979323846264338328


#include "random.c"

int main_random(int c, char *v[])
{
	float s = atof(pick_option(&c, &v, "s", "1"));
	char *offset_string = pick_option(&c, &v, "o", "");
	fprintf(stderr, "c=%d\n", c);
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s dimension distribution npoints"
		//                          0 1         2            3
		" [-s param] [-c center]\n", *v);
		return 1;
	}
	int d = atoi(v[1]);
	int n = atoi(v[3]);
	//fprintf(stderr, "generating %d %d-dimensional points from \"%s\"\n",
	//		n, d, v[2]);
	float *x = xmalloc(d*n*sizeof*x);
	double offset[d];
	for (int i = 0; i < d; i++)
		offset[i] = 0;
	if (*offset_string)
		read_n_doubles_from_string(offset, offset_string, d);

	switch (v[2][0]) {
	case 'g':
		for (int i = 0; i < n*d; i++)
			x[i] = s * random_normal();
		break;
	case 'u':
		for (int i = 0; i < n*d; i++)
			x[i] = s * (random_uniform()-0.5);
		break;
	case 'c':
		if (d != 2) fail("%d-dimensional cauchy not implemented\n", d);
		for (int i = 0; i < n; i++)
		{
			double θ = 2 * π * random_uniform();
			double ρ = s * random_cauchy();
			x[2*i+0] = ρ * cos(θ);
			x[2*i+1] = ρ * sin(θ);
		}
		break;
	default:
		return fprintf(stderr,"unrecognized distribution \"%s\"",v[2]);
	}

	for (int i = 0; i < n; i++)
	for (int j = 0; j < d; j++)
		x[i*d+j] += offset[j];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < d; j++)
			printf("%g%c", x[i*d+j], (j==d-1)?'\n':' ');

	return 0;
}



// CLI utility to access some point processing programs
int main_points(int c, char *v[])
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "random")) return main_random(c-1, v+1);
//	else if (0 == strcmp(v[1], "config")) return main_config(c-1, v+1);
//	else if (0 == strcmp(v[1], "stats")) return main_stats(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s [random|config|stats] "
			       "params... \n", *v);
	       return 1;
	}
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_points(c, v); }
#endif
