#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fail.c"
#include "xmalloc.c"
#include "colorcoords.c"
#include "drawsegment.c"
#include "marching_squares.c"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

static float get_max_length_among_small_enough(float *x, int n, float top)
{
	float r = -1;
	FORI(n) {
		float *v = x + 2*i;
		if (!isfinite(v[0]) || !isfinite(v[1]) ||
				fabs(v[0]) > top || fabs(v[1]) > top)
			continue;
		float nv = hypot(v[0], v[1]);
		if (nv > r)
			r = nv;
	}
	fprintf(stderr, "computed scale = %g\n", r);
	return r;
}

static void colorflow_flat(uint8_t *view, float *flow, int w, int h, float m)
{
	float (*x)[w][2] = (void*)flow;
	uint8_t (*y)[w][3] = (void*)view;
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		if (r > 1e8 || !isfinite(r)) {
			FORL(3) y[j][i][l] = 255; continue;
		}
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		double hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
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
	uint8_t (*i)[hack_width][3] = ii;
	i[y][x][0] = 0;
	i[y][x][1] = 0;
	i[y][x][2] = 0;
}

static void overlay_level_line_in_black(uint8_t *y,
		float *x, int w, int h, float t)
{
	int ns;
	float (*s)[2][2] = marching_squares_whole_image_float(&ns, x, w, h, t);
	black_pixel(w, h, NULL);
	for (int i = 0; i < ns; i++) {
		traverse_segment(s[i][0][0], s[i][0][1],
				s[i][1][0], s[i][1][1],
				black_pixel, y);
	}
	free(s);
}

static void overlines(uint8_t *v, float *f, int w, int h, float m)
{
	assert(m > 0);
	float *scalar = xmalloc(w * h * sizeof*scalar);
	FORI(w*h)
		scalar[i] = hypot(f[2*i], f[2*i+1]);
	FORI(51)
		overlay_level_line_in_black(v, scalar, w, h, m*i);
	free(scalar);
}

static void colorflow_ipol(uint8_t *v, float *f, int w, int h,
		float *p, int n)
{
	if (n != 1)
		fail("ipol colorization expects one single parameter");
	float satscale = p[0] ? fabs(p[0]) :
		get_max_length_among_small_enough(f, w*h, 1e5);
	colorflow_flat(v, f, w, h, satscale);
	if (p[0] < 0)
		overlines(v, f, w, h, -p[0]);
}

int middlebury_ncols = 0;
#define MIDDLEBURY_MAXCOLS 60
int middlebury_colorwheel[MIDDLEBURY_MAXCOLS][3];


void middlebury_setcols(int r, int g, int b, int k)
{
    middlebury_colorwheel[k][0] = r;
    middlebury_colorwheel[k][1] = g;
    middlebury_colorwheel[k][2] = b;
}

void middlebury_makecolorwheel(void)
{
    // relative lengths of color transitions:
    // these are chosen based on perceptual similarity
    // (e.g. one can distinguish more shades between red and yellow 
    //  than between yellow and green)
    int RY = 15;
    int YG = 6;
    int GC = 4;
    int CB = 11;
    int BM = 13;
    int MR = 6;
    middlebury_ncols = RY + YG + GC + CB + BM + MR;
    //printf("ncols = %d\n", ncols);
    if (middlebury_ncols > MIDDLEBURY_MAXCOLS)
	exit(1);
    int i;
    int k = 0;
    for (i = 0; i < RY; i++) middlebury_setcols(255, 255*i/RY, 0, k++);
    for (i = 0; i < YG; i++) middlebury_setcols(255-255*i/YG, 255, 0, k++);
    for (i = 0; i < GC; i++) middlebury_setcols(0, 255, 255*i/GC, k++);
    for (i = 0; i < CB; i++) middlebury_setcols(0, 255-255*i/CB, 255, k++);
    for (i = 0; i < BM; i++) middlebury_setcols(255*i/BM, 0, 255, k++);
    for (i = 0; i < MR; i++) middlebury_setcols(255, 0, 255-255*i/MR, k++);
}

void middlebury_computeColor(float fx, float fy, unsigned char *pix)
{
    if (middlebury_ncols == 0)
	middlebury_makecolorwheel();

    float rad = sqrt(fx * fx + fy * fy);
    float a = atan2(-fy, -fx) / M_PI;
    float fk = (a + 1.0) / 2.0 * (middlebury_ncols-1);
    int k0 = (int)fk;
    int k1 = (k0 + 1) % middlebury_ncols;
    float f = fk - k0;
    //f = 0; // uncomment to see original color wheel
    for (int b = 0; b < 3; b++) {
	float col0 = middlebury_colorwheel[k0][b] / 255.0;
	float col1 = middlebury_colorwheel[k1][b] / 255.0;
	float col = (1 - f) * col0 + f * col1;
	if (rad <= 1)
	    col = 1 - rad * (1 - col); // increase saturation with radius
	else
	    col *= .75; // out of range
	pix[2 - b] = (int)(255.0 * col);
    }
}

static bool middlebury_toolarge(float *v)
{
	return (!isfinite(v[0])) || (!isfinite(v[1])) ||
		fabs(v[0]) > 1e5 || fabs(v[1]) > 1e5;
}


static void colorflow_middlebury(uint8_t *view, float *flow, int w, int h,
		float *p, int n)
{
	if (n != 1 && n != 0)
		fail("middlebury colorization expects one single parameter");
	float (*x)[w][2] = (void*)flow;
	uint8_t (*y)[w][3] = (void*)view;
	float range = 1;
	if (n) range = p[0] ? fabs(p[0]) :
		get_max_length_among_small_enough(flow, w*h, 1e5);
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		unsigned char pix[3] = {0, 0, 0};
		if (!middlebury_toolarge(v))
			middlebury_computeColor(-v[1]/range, -v[0]/range, pix);
		FORL(3)
			y[j][i][l] = pix[l];

	}
}

static void colorflow_boldt(uint8_t *view, float *flow, int w, int h,
		float *p, int n)
{
	if (n != 1)
		fail("Homann colorization expects one single parameter");
	float (*x)[w][2] = (void*)flow;
	uint8_t (*y)[w][3] = (void*)view;
	float scale = p[0] ? fabs(p[0]) :
		get_max_length_among_small_enough(flow, w*h, 1e11)/2;
	FORJ(h) FORI(w) {
		float *vbad = x[j][i];
		float v[2] = {vbad[0]/scale, vbad[1]/scale};
		double argv = fmod((atan2(v[1], -v[0])+M_PI)/(2*M_PI),1);
		float lnv = log(1+hypot(v[0], v[1]));
		double hsv[3] = {360*argv, 1/(1+0.3*lnv), 1-1/(1.1+5*lnv)};
		double rgb[3];
		hsv_to_rgb_doubles(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];
	}
}

static void colorflow_heidelberg(uint8_t *view, float *flow, int w, int h,
		float *p, int n)
{
}

static void colorflow_depth(uint8_t *view, float *flow, int w, int h,
		float *p, int n)
{
	if (n != 1)
		fail("Depth colorization expects one single parameter");
	float (*x)[w][2] = (void*)flow;
	uint8_t (*y)[w][3] = (void*)view;
	float scale = p[0] ? fabs(p[0]) :
		get_max_length_among_small_enough(flow, w*h, 1e11)/2;
	FORJ(h) FORI(w) {
		float *vbad = x[j][i];
		float v[2] = {vbad[0]/scale, vbad[1]/scale};
		float disp = hypot(v[0]+300,v[1])-hypot(300,0);
		if (disp < 0) disp = 0;
		if (disp > 255) disp = 255;
		FORL(3)
			y[j][i][l] = disp;
	}
}

static void colorflow_shadows(uint8_t *v, float *f, int w, int h,
		float *p, int n)
{
}


static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

void colorflow(uint8_t *view, float *flow, int w, int h,
		char *method_id, char *parlist_string)
{
	float pars[0x100];
	int npars = parse_floats(pars, 0x100, parlist_string);
	if (false) {;
	} else if (0 == strcmp(method_id, "ipol")) {
		colorflow_ipol(view, flow, w, h, pars, npars);
	} else if (0 == strcmp(method_id, "middlebury")) {
		colorflow_middlebury(view, flow, w, h, pars, npars);
	} else if (0 == strcmp(method_id,"boldt")||!strcmp(method_id,"homann")){
		colorflow_boldt(view, flow, w, h, pars, npars);
	} else if (0 == strcmp(method_id, "heidelberg")) {
		colorflow_heidelberg(view, flow, w, h, pars, npars);
	} else if (0 == strcmp(method_id, "depth")) {
		colorflow_depth(view, flow, w, h, pars, npars);
	} else if (0 == strcmp(method_id, "shadows")) {
		colorflow_shadows(view, flow, w, h, pars, npars);
	} else
		fail("unrecognized colorization method %s", method_id);
}

#ifndef OMIT_MAIN
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s method pars [flow [view]]\n", *v);
				//          0 1      2     3     4
		return EXIT_FAILURE;
	}

	char *method_id = v[1];
	char *parlist_string = v[2];
	char *infile = c > 3 ? v[3] : "-";
	char *outfile = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *flow = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) fail("input is not a vector field");

	uint8_t *view = xmalloc(w * h * 3);;

	colorflow(view, flow, w, h, method_id, parlist_string);

	iio_write_image_uint8_vec(outfile, view, w, h, 3);

	free(view);
	free(flow);

	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN

// TODO:
// colorflow method "parlist" [in [out]]
//
// methods:
//
// 1. ipol
// 2. middlebury
// 3. boldt/whoever
// 4. heidelberg
// 5. depth
// 6. shadows
