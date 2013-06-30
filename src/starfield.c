#include <assert.h>
#include <math.h>
#include <stdio.h>

struct star {
	double x, y, z, b;
	int vflag;
};

struct camera {
	int w, h;
	double z, f;
	double brightness_scale;
	double brightness_min;
	double cauchy_blur;
};

// display a star in a centered camera pointing to the z direction
static int star_projection(double *oi, double *oj, double *ob,
		struct camera *c, struct star *s)
{
	double d = s->z - c->z;
	//fprintf(stderr, "%g = %g - %g\n", d, s->z, c->z);
	if (d <= 0)
		return 0;
	double F = c->f / d;
	*oi = F * s->x + c->w/2.0;
	*oj = F * s->y + c->h/2.0;
	*ob = c->brightness_scale * s->b / (d*d);
	if (*oi < 0 || *oj < 0 || *oi >= c->w || *oj >= c->h)
		return 0;
	if (*ob < c->brightness_min)
		return -1;
	return 1;
}

#include "random.c"

static double random_brightness(void)
{
	return random_pareto();
}

static double random_uniform_interval(double a, double b)
{
	return a + (b-a)*random_uniform();
}

// this function performs the inverse computation of "star projection"
static void put_visible_star(struct star *s, struct camera *c)
{
	s->b = random_brightness();
	double dtop = sqrt(s->b * c->brightness_scale / c->brightness_min);
	//fprintf(stderr, "dtop=%g\n", dtop);
	double d = random_uniform_interval(dtop/100, dtop);
	double i = random_uniform_interval(0, c->w - 1);
	double j = random_uniform_interval(0, c->h - 1);
	double F = c->f / d;
	s->z = c->z + d;
	s->x = (i - c->w/2.0)/F;
	s->y = (j - c->h/2.0)/F;
}

// this function performs the inverse computation of "star projection"
static void put_star_behind(struct star *s, struct camera *c)
{
	s->b = random_brightness();
	double dtop = sqrt(s->b * c->brightness_scale / c->brightness_min);
	double d = dtop;//random_uniform_interval(dtop/1000, dtop);
	double i = random_uniform_interval(0, c->w - 1);
	double j = random_uniform_interval(0, c->h - 1);
	double F = c->f / d;
	s->z = c->z + d;
	s->x = (i - c->w/2.0)/F;
	s->y = (j - c->h/2.0)/F;
}

#define MAIN_STARFIELD
#ifdef MAIN_STARFIELD
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"
#include "xmalloc.c"
int main(int argc, char **argv)
{
	if (argc != 6) {
		fprintf(stderr, "usage:\n\t%s w h nstars nframes dz\n", *argv);
		//                          0 1 2 3      4       5
	}
	int w = atoi(argv[1]);
	int h = atoi(argv[2]);
	int nstars = atoi(argv[3]);
	int nframes = atoi(argv[4]);
	double dz = atof(argv[5]);

	struct camera c[1];
	c->w = w;
	c->h = h;
	c->z = 0;
	c->f = 1;
	c->brightness_scale = 1;
	c->brightness_min = 1;
	c->cauchy_blur = 1;

	struct star *s = xmalloc(nstars*sizeof*s);
	for (int i = 0; i < nstars; i++)
		put_visible_star(s + i, c);

	float *x = xmalloc(w*h*sizeof*x);
	int fid = 0, frameoff = 200;
	for (int frame = 0; frame < nframes+frameoff; frame++)
	{
		memset(x, 0, w*h*sizeof*x);
		int cx = 0;
		for (int i = 0; i < nstars; i++)
		{
			double oi, oj, ob;
			int p = star_projection(&oi, &oj, &ob, c, s + i);
			if (p > 0) {
				int ii = oi;
				int jj = oj;
				assert(ii >= 0); assert(ii < w);
				assert(jj >= 0); assert(jj < h);
				float *g = x + jj*w + ii;
				*g = fmax(*g, ob);
				//x[jj*w+ii] = ob;
				cx += 1;
			} else {
				put_star_behind(s + i, c);
			}
		}
		if (frame > frameoff) {
			fprintf(stderr, "%d: %d/%d stars\n", fid,cx,nstars);

			char fnam[FILENAME_MAX];
			snprintf(fnam, FILENAME_MAX, "/tmp/xtarfield_%04d.tiff", fid);
			iio_save_image_float(fnam, x, w, h);
			fid += 1;
		}

		c->z += dz;
	}

	free(x);
	return 0;

}
#endif//MAIN_STARFIELD
