#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"
#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


static void setpixel(float *x, int w, int h, int pd, int i, int j, float *v)
{
	if (i >= 0 && j >= 0 && i < w && j < h)
		for (int l = 0; l < pd; l++)
			x[(j*w+i)*pd+l] = v[l];
}


static void overlay_one_inplace(
		float *ox, float *oy, float *of, // background images and flow
		int w, int h, int pd,            // dimensions of background
		float *ppx, int pw, int ph,      // overlaid image and its size
		float posx, float posy,          // overlay position
		float zoom, float angle,         // overlay transformation
		float dx, float dy               // overlay displacement
		)
{
	float (*px)[pw][pd] = (void*)ppx;
	float dxy[2] = {dx, dy};
	float sina = sin(angle*M_PI/180.0);
	float cosa = cos(angle*M_PI/180.0);
	for (int j = 0; j < ph; j++)
	for (int i = 0; i < pw; i++)
	{
		float cij[2] = {i - pw/2.0, j - ph/2.0};
		float ai = pw/2.0 + zoom * ( cosa * cij[0] + sina * cij[1]);
		float aj = ph/2.0 + zoom * (-sina * cij[0] + cosa * cij[1]);
		ai = round(ai);
		aj = round(aj);
		float adxy[2] = {dx + ai - i, dy + aj - j};
		// TODO: correct re-sampling (!)

		setpixel(ox,w,h,pd, posx + i      , posy + j      , px[j][i]);
		setpixel(oy,w,h,pd, posx + dx + ai, posy + dy + aj, px[j][i]);
		setpixel(of,w,h,2 , posx + i      , posy + j      , adxy);
	}
}


static void overflow_f(FILE *fs,
		float *ox, float *oy, float *of,
		float *x, float *y, float *f,
		int w, int h, int pd)
{
	// initialize output
	for (int i = 0; i < w*h*pd; i++) ox[i] = x[i];
	for (int i = 0; i < w*h*pd; i++) oy[i] = y[i];
	for (int i = 0; i < w*h* 2; i++) of[i] = f[i];

	// iterate over specified filenames
	char fname[FILENAME_MAX];
	double posx, posy, zoom, angle, dx, dy;
	while (7 == fscanf(fs, "%s %lf %lf %lf %lf %lf %lf",
				fname, &posx, &posy, &zoom, &angle, &dx, &dy)
	      )
	{
		printf("IMG=%s p=(%f %f) z=%f a=%f d=(%f %f)\n", fname,
				posx, posy, zoom, angle, dx, dy);
		int pw, ph, ppd;
		float *px = iio_read_image_float_vec(fname, &pw, &ph, &ppd);
		if (ppd != pd) fail("bad overlay pd = %d", ppd);
		overlay_one_inplace(ox, oy, of, w, h, pd,
				px, pw, ph,
				posx, posy, zoom, angle, dx, dy);
	}
	//fail("fins aquí");
}

int main(int c, char *v[])
{
	// treat input arguments
	if (c != 8) {
		fprintf(stderr, "usage:\n\t"
				"%s spec.txt bg_a bg_b bg_f o_a o_b o_f\n", *v);
				//0 1        2    3    4    5   6   7
		return 1;
	}
	char *filename_spec  = v[1];
	char *filename_bg_a  = v[2];
	char *filename_bg_b  = v[3];
	char *filename_bg_f  = v[4];
	char *filename_out_a = v[5];
	char *filename_out_b = v[6];
	char *filename_out_f = v[7];

	// read input images
	int w[3], h[3], pd[3];
	float *x = iio_read_image_float_vec(filename_bg_a, w+0, h+0, pd+0);
	float *y = iio_read_image_float_vec(filename_bg_b, w+1, h+1, pd+1);
	float *f = iio_read_image_float_vec(filename_bg_f, w+2, h+2, pd+2);
	if (pd[0] != pd[1] || pd[2] != 2)
		fail("bad pd sequence %d %d %d\n", pd[0], pd[1], pd[2]);

	// allocate output images
	float *ox = xmalloc(w[0] * h[0] * pd[0] * sizeof*ox);
	float *oy = xmalloc(w[0] * h[0] * pd[0] * sizeof*ox);
	float *of = xmalloc(w[0] * h[0] * 2 * sizeof*ox);

	// run thing
	FILE *fs = xfopen(filename_spec, "r");
	overflow_f(fs, ox, oy, of, x, y, f, *w, *h, *pd);
	xfclose(fs);

	// save output
	iio_write_image_float_vec(filename_out_a, ox, *w, *h, *pd);
	iio_write_image_float_vec(filename_out_b, oy, *w, *h, *pd);
	iio_write_image_float_vec(filename_out_f, of, *w, *h, 2);

	// cleanup and exit
	free(x); free(y); free(f);
	free(ox); free(oy); free(of);
	return 0;
}
