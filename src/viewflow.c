#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)


#include "smapa.h"
#include "fail.c"
#include "drawsegment.c"
#include "colorcoordsf.c"
#include "marching_squares.c"

static void viewflow_pd(uint8_t (**y)[3], float (**x)[2], int w, int h, float m)
{
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		float hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
		hsv_to_rgb_floats(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];

	}
}

static void viewflow_flat(uint8_t *py, float *px, int w, int h, float m)
{
	float (*x)[w][2] = (void*)px;
	uint8_t (*y)[w][3] = (void*)py;
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		if (r > 1e8 || !isfinite(r)) {
			FORL(3) y[j][i][l] = 255;
			continue;
		}
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		float hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
		hsv_to_rgb_floats(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];

	}
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


SMART_PARAMETER_SILENT(MRANGE,0)

static bool middlebury_toolarge(float *v)
{
	return (!isfinite(v[0])) || (!isfinite(v[1])) ||
		fabs(v[0]) > 1e5 || fabs(v[1]) > 1e5;
}

static void viewflow_middlebury(uint8_t *py, float *px, int w, int h)
{
	float (*x)[w][2] = (void*)px;
	uint8_t (*y)[w][3] = (void*)py;
	float range = MRANGE();
	if (!isfinite(range)) {
		range = -1;
		FORJ(h) FORI(w) {
			float *v = x[j][i];
			if (middlebury_toolarge(v))
				continue;
			float nv = hypot(v[0], v[1]);
			if (nv > range)
				range = nv;
		}
	}
	fprintf(stderr, "range = %g\n", range);
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		unsigned char pix[3] = {0, 0, 0};
		if (!middlebury_toolarge(v))
			middlebury_computeColor(-v[1]/range, -v[0]/range, pix);
		FORL(3)
			y[j][i][l] = pix[l];

	}
}

static float pick_scale(float (*x)[2], int n)
{
	float range = -1;
	FORI(n) {
		if (middlebury_toolarge(x[i]))
			continue;
		float nv = hypot(x[i][0], x[i][1]);
		if (nv > range)
			range = nv;
	}
	return range;
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

SMART_PARAMETER_SILENT(NOVERLINES,51)

static void *matrix_build(int w, int h, size_t n)
{
	size_t p = sizeof(void *);
	char *r = xmalloc(h*p + w*h*n);
	for (int i = 0; i < h; i++)
		*(void **)(r + i*p) = r + h*p + i*w*n;
	return r;
}



static void overlines(uint8_t (**y)[3], float (**x)[2], int w, int h, float s)
{
	assert(s > 0);
	float **scalar = matrix_build(w, h, sizeof**scalar);
	FORJ(h) FORI(w)
		scalar[j][i] = hypot(x[j][i][0], x[j][i][1]);

	int nl = NOVERLINES();
	FORI(nl)
		overlay_level_line_in_black(y, scalar, w, h, s*i);
}

static char *help_string_name     = "viewflow";
static char *help_string_version  = "viewflow 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "represent a vector field using a color code";
static char *help_string_usage    = "usage:\n\t"
"SATSCALE [in.flo [out.png]]";
static char *help_string_long     =
"Viewflow displays a two-dimensional vector field using a color.\n"
"\ni.e. "
"The input is an image with two-dimensional pixels (i.e. a vector field).\n"
"The output is an RGB image of the same size, where each vector is\n"
"represented by a different color, according to a two-dimensional palette.\n"
"For most palettes, the hue indicates the direction and the intensity\n"
"indicates the magnitude of the vector.\n"
"Vectors larger than SATSCALE are \"saturated\" to the same intensity.\n"
"\n"
"Usage: viewflow SATSCALE in.flo out.png\n"
"   or: viewflow SATSCALE in.flo > out.npy\n"
"   or: cat in.flo | viewflow SATSCALE > out.npy\n"
"\n"
"Environment:\n"
" MRANGE\tMaximum range for Middlebury palette (default 0)\n"
" NOVERLINES\tTotal number of level lines to draw (default 50)\n"
"Options:\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
"\n"
"Examples:\n"
" plambda a.png x,g | viewflow 50  View gradient of grayscale image\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above
int main_viewflow(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);

	if (c != 4 && c != 3 && c != 2) {
		fprintf(stderr, "usage:\n\t%s satscale flow view\n", *v);
				//          0    1       2     3
		return EXIT_FAILURE;
	}

	char *infile = c > 2 ? v[2] : "-";
	char *outfile = c > 3 ? v[3] : "-";

	int w, h, pd;
	void *data = iio_read_image_float_matrix_vec(infile, &w, &h, &pd);
	float (**flow)[pd] = data;
	if (pd != 2) fail("input is not a vector field");

	uint8_t (**view)[3] = matrix_build(w, h, sizeof**view);

	float satscale = atof(v[1]);
	if (0 == satscale) {
		satscale = pick_scale(flow[0], w*h);
		//fprintf(stderr, "computed scale = %g\n", satscale);
	}

	//viewflow_pd(view, flow, w, h, fabs(satscale));
	if (isfinite(satscale)) {
		viewflow_flat(view[0][0], flow[0][0], w, h, fabs(satscale));
		if (satscale < 0)
			overlines(view, flow, w, h, -satscale);
	} else
		viewflow_middlebury(view[0][0], flow[0][0], w, h);

	iio_write_image_uint8_vec(outfile, view[0][0], w, h, 3);

	free(view);
	free(flow);

	return EXIT_SUCCESS;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char *v[]) { return main_viewflow(c, v); }
#endif//

// TODO:
// viewflow method "parlist" [in [out]]
//
// methods:
//
// 1. ipol
// 2. middlebury
// 3. boldt/whoever
// 4. heidelberg
// 5. depth
// 6. shadows
