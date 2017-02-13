#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#ifdef IIO_MAX_DIMENSION
# define MAX_PIXELDIM IIO_MAX_DIMENSION
#else
# define MAX_PIXELDIM 20
#endif

struct vecstats {
	int pd;
	float min[MAX_PIXELDIM];
	float max[MAX_PIXELDIM];
	float avg[MAX_PIXELDIM];
};

static float euclidean_norm(float *x, int n)
{
	assert(n > 0);
	return n-1 ? hypot(*x, euclidean_norm(x+1, n-1)) : fabs(*x);
}

static void get_vecstats(struct vecstats *s, float *x, int n, int pd)
{
	s->pd = pd;
	float min = INFINITY;
	float max = -INFINITY;
	FORL(pd) s->avg[l] = 0;

	FORI(n) {
		float t = euclidean_norm(x + i*pd, pd);
		if (t < min) {
			min = t;
			FORL(pd) s->min[l] = x[i*pd + l];
		}
		if (t > max) {
			max = t;
			FORL(pd) s->max[l] = x[i*pd + l];
		}
		FORL(pd)
			s->avg[l] += x[i*pd + l] / n;
	}
}

static void print_vec(float *x, int n)
{
	assert(n > 0);
	printf("(%g", x[0]);
	for (int i = 1; i < n; i++)
		printf(",%g", x[i]);
	printf(")");
}

static void print_vecstats(struct vecstats *s)
{
	if (s->pd == 1)
		printf("min=%g\navg=%g\nmax=%g\n",s->min[0],s->avg[0],s->max[0]);
	else {
		printf("min="); print_vec(s->min, s->pd);
		printf("\navg="); print_vec(s->avg, s->pd);
		printf("\nmax="); print_vec(s->max, s->pd);
		printf("\n");
	}
}

int main(int c, char *v[])
{
	if (c != 1 && c != 2) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		return EXIT_FAILURE;
	}
	char *filename = c > 1 ? v[1] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename, &w, &h, &pd);

	printf("%dx%d", w, h);
	if (pd > 1) printf("x%d", pd);
	//printf("  ");
	printf("\n");
	struct vecstats s[1];
	get_vecstats(s, x, w*h, pd);
	print_vecstats(s);

	return EXIT_SUCCESS;
}
