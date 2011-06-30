#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"
#include "marching_squares.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)



#include "fragments.c"

static void viewflow_pd(uint8_t (**y)[3], float (**x)[2], int w, int h, float m)
{
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		double hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
		void hsv_to_rgb_doubles(double*, double*);
		hsv_to_rgb_doubles(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];

	}
}

static void viewflow_flat(uint8_t *py, float *px, int w, int h, int m)
{
	float (*x)[w][2] = (void*)px;
	uint8_t (*y)[w][3] = (void*)py;
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		double hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
		void hsv_to_rgb_doubles(double*, double*);
		hsv_to_rgb_doubles(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];

	}
}

static void black_pixel(int x, int y, void *ii)
{
	static int hack_width = 0;
	static int hack_height = 0;
	if (!ii) {
		hack_width = x;
		hack_height = y;
		return;
	}
	if (x < 0 || y < 0 || x > hack_width || y > hack_height) return;
	uint8_t (**i)[3] = ii;
	i[y][x][0] = 0;
	i[y][x][1] = 0;
	i[y][x][2] = 0;
}

static void overlay_level_line_in_black(uint8_t (**y)[3],
		float **x, int w, int h, float t)
{
	int ns;
	float (*s)[2][2] = marching_squares_whole_image_float(&ns, *x, w, h, t);
	black_pixel(w, h, NULL);
	for (int i = 0; i < ns; i++) {
		traverse_segment(s[i][0][0], s[i][0][1],
				s[i][1][0], s[i][1][1],
				black_pixel, y);
	}
	free(s);
}

static void overlines(uint8_t (**y)[3], float (**x)[2], int w, int h, float s)
{
	assert(s > 0);
	float **scalar = matrix_build(w, h, sizeof**scalar);
	FORJ(h) FORI(w)
		scalar[j][i] = hypot(x[j][i][0], x[j][i][1]);

	FORI(11)
		overlay_level_line_in_black(y, scalar, w, h, s*i);
}

int main_viewflow(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s satscale flow view\n", *v);
				//          0    1       2     3
		return EXIT_FAILURE;
	}

	int w, h, pd;
	void *data = iio_read_image_float_matrix_vec(v[2], &w, &h, &pd);
	float (**flow)[pd] = data;
	if (pd != 2) error("input is not a vector field");

	uint8_t (**view)[3] = matrix_build(w, h, sizeof**view);

	float satscale = atof(v[1]);

	//viewflow_pd(view, flow, w, h, fabs(satscale));
	viewflow_flat(view[0][0], flow[0][0], w, h, fabs(satscale));
	if (satscale < 0)
		overlines(view, flow, w, h, -satscale);

	iio_save_image_uint8_vec(v[3], view[0][0], w, h, 3);

	free(view);
	free(flow);

	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_viewflow(c, v);
}
#endif//OMIT_MAIN
