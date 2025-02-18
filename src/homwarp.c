
// y = H(x)
static void apply_homography(double y[2], double H[9], double x[2])
{
	double X = H[0]*x[0] + H[1]*x[1] + H[2];
	double Y = H[3]*x[0] + H[4]*x[1] + H[5];
	double Z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = X / Z;
	y[1] = Y / Z;
}

static double invert_homography(double o[9], double i[9])
{
    double det = i[0]*i[4]*i[8] + i[2]*i[3]*i[7] + i[1]*i[5]*i[6]
               - i[2]*i[4]*i[6] - i[1]*i[3]*i[8] - i[0]*i[5]*i[7];
    o[0] = (i[4]*i[8] - i[5]*i[7]) / det;
    o[1] = (i[2]*i[7] - i[1]*i[8]) / det;
    o[2] = (i[1]*i[5] - i[2]*i[4]) / det;
    o[3] = (i[5]*i[6] - i[3]*i[8]) / det;
    o[4] = (i[0]*i[8] - i[2]*i[6]) / det;
    o[5] = (i[2]*i[3] - i[0]*i[5]) / det;
    o[6] = (i[3]*i[7] - i[4]*i[6]) / det;
    o[7] = (i[1]*i[6] - i[0]*i[7]) / det;
    o[8] = (i[0]*i[4] - i[1]*i[3]) / det;
    return det;
}

static void matrix_product(double o[9], double x[9], double y[9])
{
	// x0 x1 x2   y0 y1 y2
	// x3 x4 x5   y3 y4 y5
	// x6 x7 x8   y6 y7 y8

	double t[9];

	t[0] = x[0]*y[0] + x[1]*y[3] + x[2]*y[6];
	t[1] = x[0]*y[1] + x[1]*y[4] + x[2]*y[7];
	t[2] = x[0]*y[2] + x[1]*y[5] + x[2]*y[8];

	t[3] = x[3]*y[0] + x[4]*y[3] + x[5]*y[6];
	t[4] = x[3]*y[1] + x[4]*y[4] + x[5]*y[7];
	t[5] = x[3]*y[2] + x[4]*y[5] + x[5]*y[8];

	t[6] = x[6]*y[0] + x[7]*y[3] + x[8]*y[6];
	t[7] = x[6]*y[1] + x[7]*y[4] + x[8]*y[7];
	t[8] = x[6]*y[2] + x[7]*y[5] + x[8]*y[8];

	for (int i = 0; i < 9; i++)
		o[i] = t[i];
}


#include <math.h>
//#include "extrapolators.c"
//#include "bilinear_interpolation.c"
//#include "marching_interpolation.c"
//#include "bicubic_gray.c"
#include "spline.c"

// SECTION 4. Extrapolation                                                 {{{1

// A "extrapolator" evaluates an image at an arbitrary integral position.
// When the position is outside the image domain, the value is extrapolated
// (by periodization, or by a constant value).

// type of the "extrapolator" functions
typedef float (*extrapolator_t)(float*,int,int,int,int,int,int);

// auxiliary function: compute n%p correctly, even for huge and negative numbers
//static int good_modulus(int n, int p)
//{
//	int r = n % p;
//	r = r < 0 ? r + p : r;
//	assert(r >= 0);
//	assert(r < p);
//	return r;
//}

// instance of "extrapolator_t", extrapolate by periodicity
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// instance of "extrapolator_t", extrapolate by a constant value
static float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return value;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

//// obtain an extrapolator function from the current configuration
//static extrapolator_t obtain_extrapolator(struct viewer_state *e)
//{
//	if (e->tile_plane)
//		return getsample_per;
//	float top = 255; // white
//	getsample_cons(&top, 0, 0, 0, 0, 0, 0);
//	return getsample_cons;
//}



// SECTION 5. Local Interpolation                                           {{{1
//
// An "interpolator" evaluates an image at an arbitrary floating-point
// position.  When the position is not integral, the value is a combination of
// neighboring values.  Notice that when the position is outside the image
// domain, an extrapolator is also needed.

// type of the "interpolator" functions
typedef float (*interpolator_t)(float*,int,int,int,float,float,int,
		extrapolator_t);

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// instance of "interpolator_t", for bilinear interpolation
static float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}

// auxiliary code for "linear" interpolation
#include "marching_interpolation.c"

// instance of "interpolator_t", for linear (marching) interpolation
static float linear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return marchi(a, c, b, d, p-ip, q-iq);
}

// instance of "interpolator_t" for nearest neighbor interpolation
static float nearest_neighbor_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = round(p);
	int iq = round(q);
	return pix(x, w, h, pd, ip, iq, l);
}


// one-dimensional cubic interpolation of four data points ("Keys")
static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

// two-dimensional separable cubic interpolation, on a 4x4 grid
static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

// instance of "interpolator_t" for bicubic interpolation
static float bicubic_interpolation_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}

//// obtain an interpolation operator from the current configuration
//static interpolator_t obtain_interpolator(struct viewer_state *e)
//{
//	if (e->dragged_point >= 0 || e->dragging_background)
//		return nearest_neighbor_at;
//	if (e->interpolation_order == 1) return linear_interpolation_at;
//	if (e->interpolation_order == 2) return bilinear_interpolation_at;
//	if (e->interpolation_order == 3) return bicubic_interpolation_at;
//	return nearest_neighbor_at;
//}

#include "parsenumbers.c"
static void read_homography_from_string(double H[9], char *s)
{
	H[0] = 1; H[1] = 0; H[2] = 0;
	H[3] = 0; H[4] = 1; H[5] = 0;
	H[6] = 0; H[7] = 0; H[8] = 1;
	double x[2];
	if (1 == sscanf(s, "r%lg", x)) {
		double a = M_PI**x/180;
		fprintf(stderr, "rotation of angle %g (%g)\n", *x, a);
		H[0] = H[4] = cos(a);
		H[1] = -sin(a);
		H[3] = sin(a);
		return;
	}
	read_n_doubles_from_string(H, s, 9);
}


// SECTION 6. Main warping function                                        {{{1


int homwarp(float *X, int W, int H, double M[9], float *x,
		int w, int h, int o, int bd)
{
	interpolator_t u = bicubic_interpolation_at;
	if (o == 0) u = nearest_neighbor_at;
	if (o == 1) u = linear_interpolation_at;
	if (o == 2) u = bilinear_interpolation_at;
	//if (o ==-2) u = quilez3_interpolation_at;
	//if (o ==-4) u = quilez5_interpolation_at;

	extrapolator_t e = bd ? getsample_per : getsample_cons;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, M, p);
		X[j*W+i] = u(x, w, h, 1, p[0], p[1], 0, e);
	}

	return 0;
}

int shomwarp(float *X, int W, int H, double M[9], float *x,
		int w, int h, int o, int bd)
{
	// if low order-interpolation, evaluate right away
	if (o == 0 || o == 1 || o == 2 || o == -2 || o == -3)
		return homwarp(X, W, H, M, x, w, h, o, bd);

	// otherwise, pre-filtering is required
	bool r = prepare_spline(x, w, h, 1, o);
	if (!r) return 2;

	// warp the points
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, M, p);
		p[0] += 0.5; // solve a mis-alignement convention
		p[1] += 0.5;
		float *out = X + (j*W + i);
		evaluate_spline_at(out, x, w, h, 1, o, p[0], p[1]);
	}

	return 0;
}

// now begins the main function of the CLI interface
static char *help_string_name     = "homwarp";
static char *help_string_version  = "homwarp 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "warp an image by an homography";
static char *help_string_usage    = "usage:\n\t"
"homwarp [-i] [-o ORDER] HOMOGRAPHY [W H [in [out]]]";
static char *help_string_long     =
"Homwarp resamples an image according to an homography.\n"
"\n"
"A new image is created whose pixel values are sampled at the positions\n"
"specified by an homographic transform.  No fancy filtering is performed,\n"
"thus this program is not appropriate for zooming-out images as it will\n"
"introduce a lot of aliasing.\n"
"\n"
"The only required parameter is HOMOGRAPHY, which is a string that specifies\n"
"the homographic transform (either a single string with 8 numbers or a\n"
"filename that contains them).  Optionally, the user can also specify the\n"
"desired width and height of the output image (by default, the same as the\n"
"input image), and the input/output files (by default, stdin/stdout).\n"
"\n"
"Usage: homwarp HOMOGRAPHY < in > out\n"
"   or: homwarp HOMOGRAPHY WIDTH HEIGHT < in > out\n"
"   or: homwarp HOMOGRAPHY WIDTH HEIGTH in > out\n"
"   or: homwarp HOMOGRAPHY WIDTH HEIGTH in out\n"
"\n"
"Options:\n"
" -i        use the inverse homography\n"
" -p        use periodic boundary conditions\n"
" -o ORDER  choose a different interpolation method (default = -3)\n"
" -h        display short help message\n"
" --help    display longer help message\n"
"\n"
"Homography:\n"
" \"a b p c d q r s t\"\tread the 9 coefficients of the matrix from a string\n"
" homography.txt\tread the 9 coefficients of the matrix from a file\n"
"\n"
"Order:\n"
" 0\tnearest neighbor interpolation\n"
" 1\tpiecewise linear interpolation\n"
" 2\tbilinear interpolation\n"
" -3\tbicubic interpolation\n"
//" -2\tcubic faded bilinear interpolation (a la Inigo Quilez)\n"
//" -4\tquintic faded bilinear interpolation (a la Inigo Quilez)\n"
" 3,5,7\tspline of order 3,5 or 7\n"
"\n"
"Examples:\n"
" homwarp \"1 0 X 0 1 Y 0 0 1\" W H < in > out\tcrop at (X,Y) of size WxH\n"
" homwarp \"-1 0 640 0 1 0 0 0 1\" < in > out\thoriz. flip of a 640x480 image\n"
" homwarp \"[-1,0,640;0,1,0;0,0,1];\" < in > out\tnon-numeric chars. are ignored\n"
" homwarp homography.txt < in > out          \tread homography from file\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
#include "help_stuff.c"
#include "pickopt.c"
int main_homwarp(int c, char *v[])
{
	if (c == 2) if_help_is_requested_print_it_and_exit_the_program(v[1]);
	bool do_invert = pick_option(&c, &v, "i", NULL);
	bool do_center = pick_option(&c, &v, "c", NULL);
	int ord = atoi(pick_option(&c, &v, "o", "-3"));
	int bd = !!pick_option(&c, &v, "p", NULL);
	if (c != 2 && c != 4 && c != 5 && c != 6)
		return fprintf(stderr, "usage:\n\t"
		"%s [-i] [-o {0|1|2|-3|3|5|7}] hom [w h [in [out]]]\n"
		//0                            1   2 3  4   5
		"\t-i\tinvert input homography\n"
		"\t-o\tchose interpolation order (default -3 = bicubic)\n"
		, *v);
	double H_direct[9], H_inv[9];
	read_homography_from_string(H_direct, v[1]);
	invert_homography(H_inv, H_direct);
	double *H = do_invert ? H_inv : H_direct;
	int ow = c > 2 ? atoi(v[2]) : -1;
	int oh = c > 3 ? atoi(v[3]) : -1;
	char *filename_in  = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	if (ow == -1) ow = w;
	if (oh == -1) oh = h;
	float *y = xmalloc(ow * oh * pd * sizeof*y);

	if (do_center) {
		double P[9] = { 1, 0, -w/2, 0, 1, -h/2, 0, 0, 1};
		double iP[9] = { 1, 0, w/2, 0, 1, h/2, 0, 0, 1};
		matrix_product(H, H, P);
		matrix_product(H, iP, H);
	}

	int r = 0;
	for (int i = 0; i < pd; i++)
		r += shomwarp(y + i*ow*oh, ow, oh, H, x + i*w*h, w, h, ord, bd);

	iio_write_image_float_split(filename_out, y, ow, oh, pd);
	return r;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_homwarp(c, v); }
#endif
