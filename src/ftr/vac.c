// vac: visualisateur d'autocorrélation


// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h> // for getpid only

// user interface library
#include "ftr.h"

// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


#include "random.c"

#define OMIT_BLUR_MAIN
#include "blur.c"

#define OMIT_PPSMOOTH_MAIN
#include "ppsmooth.c"


// fourier stuff

// data structure to store the state of the viewer
struct viewer_state {
	// structural parameters
	int w;  // scene width
	int h;  // scene height
	int o;  // scene orientation 0=horizontal, 1=vertical

	// input texture parameters
	int d;          // random seed
	float g;        // gaussian grain of the texture
	int l;          // number of laplacian steps
	int s[2][3];    // two shifts (x,y,sign)
	float p;        // saturation parameter (percentile)

	// autocorrelation parameters
	float r;        // center radius mask
	float z;        // saturation parameter

	// ui
	struct bitmap_font font[1];
};


// function to reset and center the viewer
static void center_state(struct viewer_state *e)
{
	// input texture
	e->d = 1;
	e->g = 4;
	e->l = 1;
	e->s[0][0] = 100; e->s[0][1] =  0; e->s[0][2] = +1;
	e->s[1][0] =  50; e->s[1][1] = 87; e->s[1][2] = -1;
	e->p = 1;

	// autocorrelation view
	e->r = 20;
	e->z = 0.1;
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}


// SECTION 2. algorithms                                                    {{{1

typedef float (*getsample_operator)(float*,int,int,int,int);

// like n%p, but works for all numbers
static int gmod(int x, int m)
{
	int r = x % m;
	return r < 0 ? r + m : r;
}

// extrapolate by periodicity
inline static float getsample_periodic(float *x, int w, int h, int i, int j)
{
	i = gmod(i, w);
	j = gmod(j, h);
	return x[j*w+i];
}

static void fill_base_texture(float *x, struct viewer_state *e)
{
	int w = e->w;
	int h = e->h;

	// create gaussian noise with seed e->d
	xsrand(e->d);
	for (int i = 0; i < w * h; i++)
		x[i] = random_normal();

	// blur it with grain e->g
	float g[1] = { e->g };
	blur_2d(x, x, w, h, 1, "gaussian", g, 1);

	// apply the given shifts
	getsample_operator P = getsample_periodic;
	float *t = malloc(w * h * sizeof*t);
	for (int i = 0; i < w * h; i++)
		t[i] = x[i];
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		x[j*w + i] += e->s[0][2] * P(t,w,h, i+e->s[0][0], j+e->s[0][1]);
		x[j*w + i] += e->s[1][2] * P(t,w,h, i+e->s[1][0], j+e->s[1][1]);
	}
	free(t);
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void getpercentiles(float *m, float *M, float *x, int n, float q)
{
	float *t = xmalloc(n * sizeof*t);
	memcpy(t, x, n*sizeof*t);
	int N = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(x[i]))
			t[N++] = x[i];
	//fprintf(stderr, "getperc n=%d N=%d p=%g\n", n, N, N*100.0/n);
	qsort(t, N, sizeof*t, compare_floats);
	int a = q/100*N;
	int b = (1-q/100)*N;
	if (a < 0) a = 0;
	if (b < 0) b = 0;
	if (a >= N) a = N-1;
	if (b >= N) b = N-1;
	*m = t[a];
	*M = t[b];
	free(t);
}

static unsigned char bclamp(float x)
{
	if (x < 0) return 0;
	if (x > 255) return 255;
	return x;
}

static void display_base_texture(uint8_t *vx, float *x, struct viewer_state *e)
{
	// compute e->p percentiles
	float m, M;
	getpercentiles(&m, &M, x, e->w * e->h, e->p);

	// saturate at these percentiles
	for (int i = 0; i < e->w * e->h; i++)
		vx[i] = bclamp( 255 * (x[i] - m) / (M - m) );
}

static void compute_autocorrelation(float *y, float *x, int w, int h)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *ys = xmalloc(w*h*sizeof*c);
	fftwf_complex *fc = fftwf_xmalloc(w*h*sizeof*fc);
	for (int l = 0; l < 1; l++)
	{
		for (int i = 0; i < w*h; i++)
			c[i] = x[i];
		fft_2dfloat(fc, c, w, h);
		for (int i = 0; i < w*h; i++)
			fc[i] = cabs(fc[i]);
		ifft_2dfloat(c, fc, w, h);
		for (int i = 0; i < w*h; i++)
			ys[i] = c[i];
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			int ii = (i + w/2) % w;
			int jj = (j + h/2) % h;
			y[j*w+i] = ys[jj*w+ii] * 1;//norm;
		}
	}
	free(ys);
	free(c);
	fftwf_free(fc);
}

static
void display_autocorrelation(uint8_t *vX, float *X, struct viewer_state *e)
{
	//for (int i = 0; i < 3 * e->w * e->h; i++)
	//	vX[i] = 127;

	// compute e->z percentiles
	float m, M;
	getpercentiles(&m, &M, X, e->w * e->h, e->z);

	// saturate at these percentiles
	for (int i = 0; i < e->w * e->h; i++)
	{
		vX[3*i+0] = bclamp( 255 * (X[i] - m) / (M - m) );
		vX[3*i+1] = bclamp( 255 * (X[i] - m) / (M - m) );
		vX[3*i+2] = bclamp( 255 * (X[i] - m) / (M - m) );
	}
}




// SECTION 3. Coordinate Conversions                                        {{{1



// SECTION 4. Drawing                                                       {{{1





// Subsection 7.2. Drawing user-interface elements                          {{{2


// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas (dark blue)
	for (int i = 0 ; i < f->w * f->h; i++)
	{
		f->rgb[3*i+0] = 0;
		f->rgb[3*i+1] = 0;
		f->rgb[3*i+2] = 100;
	}

	// space for images and their visualizations
	float *x = xmalloc(e->w * e->h * sizeof*x);
	float *X = xmalloc(e->w * e->h * sizeof*X);
	uint8_t *vx = xmalloc(e->w * e->h);
	uint8_t *vX = xmalloc(3 * e->w * e->h);

	// perform all the computaitons
	fill_base_texture(x, e);
	display_base_texture(vx, x, e);
	compute_autocorrelation(X, x, e->w, e->h);
	display_autocorrelation(vX, X, e);

	// window 0: base texture
	// WxH starting at (0,0)
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
	{
		uint8_t g = vx[j*e->w + i];
		f->rgb[(j*f->w + i)*3 + 0] = g;
		f->rgb[(j*f->w + i)*3 + 1] = g;
		f->rgb[(j*f->w + i)*3 + 2] = g;
	}

	// window 1: autocorrelation
	// if o==0 : WxH starting at (w,0)
	// if o==1 : WxH starting at (0,h)
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
	{
		uint8_t *c = vX + 3*(j*e->w + i);
		int oi = e->w * (1 - e->o);
		int oj = e->h * e->o;
		f->rgb[((j+oj)*f->w + i+oi)*3 + 0] = c[0];
		f->rgb[((j+oj)*f->w + i+oi)*3 + 1] = c[1];
		f->rgb[((j+oj)*f->w + i+oi)*3 + 2] = c[2];
	}


	// free image data
	free(X); free(vX);
	free(x); free(vx);


	// hud
	//uint8_t fg[3] = {0, 0, 0};
	//uint8_t bg[3] = {127, 127, 127};
	uint8_t fg[3] = {0, 255, 0};
	//uint8_t bg[3] = {100, 100, 100};
	uint8_t bg[3] = {0, 0, 0};
	//uint8_t fg[3] = {255, 0, 255};
	//uint8_t *bg = 0;
	char buf[0x200];
	snprintf(buf, 0x200,
			"d (random seed)    = %d\n"
			"g (gaussian grain) = %g\n"
			"l (num laplacians) = %d\n"
			"s1 (shift 1)       = %3d %3d %3d\n"
			"s2 (shift 2)       = %3d %3d %3d\n"
			"p (texture saturation)    = %g\n"
			"r (autocorr. mask radius) = %g\n"
			"z (autocorr. saturation)  = %g\n",
			e->d, e->g, e->l,
			e->s[0][0], e->s[0][1], e->s[0][2],
			e->s[1][0], e->s[1][1], e->s[1][2],
			e->p, e->r, e->z);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0+0, 0+0, fg, bg, 0, e->font, buf);
}



// SECTION 5. User-Interface Actions and Events                             {{{1


//// action: scale strata frequency
//static void scale_strata_frequency(struct viewer_state *e, float f)
//{
//	e->f *= f;
//	//fprintf(stderr, "f = %g\n", e->f);
//}
//static void scale_fold_parameter(struct viewer_state *e, float f)
//{
//	e->p *= f;
//	//fprintf(stderr, "p = %g\n", e->p);
//}
//static void shift_angle_a(struct viewer_state *e, float s) { e->a += s; }
//static void shift_angle_b(struct viewer_state *e, float s) { e->b += s; }
//static void shift_angle_c(struct viewer_state *e, float s) { e->c += s; }
//static void scale_radius(struct viewer_state *e, float f) { e->R *= f; }
//static void shift_shift_s(struct viewer_state *e, float s) { e->s += s; }
//static void shift_shift_z(struct viewer_state *e, float s) { e->z0 += s; }

static void shift_random_seed(struct viewer_state *e, float s) { e->d += s; }
static void scale_gaussian_grain(struct viewer_state *e, float f) { e->g *= f; }

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	int p = getpid();
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_vac_%d_%d.png", p, c);
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(n, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}


// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (islower(k) && m&FTR_MASK_SHIFT)
		k = toupper(k);
	//fprintf(stderr, "key k=%d m=%d\n", k, m);

	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 1);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'z') center_view(f);
//	if (k == 'f') scale_strata_frequency(e, 1.3);
//	if (k == 'F') scale_strata_frequency(e, 1/1.3);
//	if (k == 'p') scale_fold_parameter(e, 1.3);
//	if (k == 'P') scale_fold_parameter(e, 1/1.3);
////	if (k == 's') shift_shift_s(e, 1);
////	if (k == 'S') shift_shift_s(e, -1);
//	if (k == 'r') scale_radius(e, 1.3);
//	if (k == 'R') scale_radius(e, 1/1.3);
//	if (k == 'a') shift_angle_a(e, 10);
//	if (k == 'A') shift_angle_a(e, -10);
//	if (k == 'b') shift_angle_b(e, 2);
//	if (k == 'B') shift_angle_b(e, -2);
//	if (k == 'c') shift_angle_c(e, 2);
//	if (k == 'C') shift_angle_c(e, -2);
	if (k == ',') action_screenshot(f);

	f->changed = 1;
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	f->changed = 1;
}

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// d, g, l, b, c, R, s hitboxes of font height
	// 0  1  2  3  4  5  6
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) shift_random_seed(e, -1);
		if (Y == 1) scale_gaussian_grain(e, 1/1.3);
//		if (Y == 1) scale_fold_parameter(e, 1/1.3);
//		if (Y == 2) shift_angle_a(e, -10);
//		if (Y == 3) shift_angle_b(e, -2);
//		if (Y == 4) shift_angle_c(e, -2);
//		if (Y == 5) scale_radius(e, 1/1.3);
//		if (Y == 6) shift_shift_s(e, -5);
//		if (Y == 7) shift_shift_z(e, -5);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) shift_random_seed(e, 1);
		if (Y == 1) scale_gaussian_grain(e, 1.3);
//		if (Y == 1) scale_fold_parameter(e, 1.3);
//		if (Y == 2) shift_angle_a(e, 10);
//		if (Y == 3) shift_angle_b(e, 2);
//		if (Y == 4) shift_angle_c(e, 2);
//		if (Y == 5) scale_radius(e, 1.3);
//		if (Y == 6) shift_shift_s(e, 5);
////		if (Y == 7) shift_shift_z(e, 5);
	}

	f->changed = 1;
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	;
}

// expose handler
static void event_expose(struct FTR *f, int b, int m, int x, int y)
{
	if (f->changed)
		paint_state(f);
}



// SECTION 6. Main Program                                                 {{{1

#include "pickopt.c"

//int main_vac_noninteractive(int c, char *v[])
//{
//	// initialize state (sets dummy default arguments)
//	struct viewer_state e[1];
//	center_state(e);
//	e->stratum = 0;
//
//	// process named arguments
//	int w = atoi(pick_option(&c, &v, "w", "360"));
//	int h = atoi(pick_option(&c, &v, "h", "720"));
//	char *filename_out = pick_option(&c, &v, "o", "-");
//	char *filename_stratum = pick_option(&c, &v, "S", "");
//	if (*filename_stratum)
//	{
//		uint8_t *iio_read_image_uint8(char*,int*,int*);
//		e->stratum = iio_read_image_uint8(filename_stratum,
//				&e->stratum_w, &e->stratum_h);
//		e->f = 1;
//	}
////	e->f = atof(pick_option(&c, &v, "f", "0.1"));
////	e->p = atof(pick_option(&c, &v, "p", "0.001"));
////	e->a = atof(pick_option(&c, &v, "a", "0"));
////	e->b = atof(pick_option(&c, &v, "b", "0"));
////	e->c = atof(pick_option(&c, &v, "c", "0"));
////	e->R = atof(pick_option(&c, &v, "R", "50"));
//
//	// fill-in output image
//	uint8_t *x = malloc(3*w*h);
//	paint_cylinder(x, w, h, e);
//
//	// save output image
//	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
//	iio_write_image_uint8_vec(filename_out, x, w, h, 3);
//
//	return 0;
//}

static char *help_string_name     = "vac";
static char *help_string_version  = "vac 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "visualize auto correlations";
static char *help_string_usage    = "usage:\n\tvac";
static char *help_string_long     =
"Vac is an interface for exploring autocorrelation images.\n"
"\n"
"A geological structure in the shape of a parallel parabolic fold\n"
"is traversed by a cylindrical borehole.  You can rotate the cylinder\n"
"and look at the intersection of the strata.  To change parameters,\n"
"use the mouse wheel over each parameter name."
"\n"
"Usage: vac\n"
"   or: vac -n [options]\n"
"\n"
"Keys:\n"
" q,ESC  quit the program\n"
" z      go back to default parameters\n"
"\n"
"Options:\n"
"\n"
" -n\tenable non-interactive mode\n"
" -o X\twrite output to file X (default stdout)\n"
" -w X\tset image width (default=800)\n"
" -h X\tset image height (default=800)\n"
" -f X\tset strata frequency parameter (default=0.1)\n"
" -p X\tset parabola curvature (default=0.001)\n"
" -a X\tset cylinder orientation, first euler angle (default=0)\n"
" -b X\tset cylinder orientation, second euler angle (default=0)\n"
" -c X\tset cylinder orientation, third euler angle (default=0)\n"
" -R X\tset cylinder radius (default=50)\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>.\n"
;

#include "help_stuff.c"
int main_vac(int c, char *v[])
{
	// if requested, print help
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

//	// if -n option, run noninteractively
//	if (pick_option(&c, &v, "n", 0))
//		return main_vac_noninteractive(c, v);

	// initialize state
	struct viewer_state e[1];
	int N = 400;
	e->w = N;
	e->h = N;
	e->o = 0;
	center_state(e);
//	e->stratum = 0;

	// extract named arguments
	char *filename_out = pick_option(&c, &v, "o", "-");
//	char *filename_stratum = pick_option(&c, &v, "S", "");
//	if (*filename_stratum)
//	{
//		uint8_t *iio_read_image_uint8(char*,int*,int*);
//		e->stratum = iio_read_image_uint8(filename_stratum,
//				&e->stratum_w, &e->stratum_h);
//		e->f = 1;
//	}
//	e->f = atof(pick_option(&c, &v, "f", "0.1"));
//	e->p = atof(pick_option(&c, &v, "p", "0.001"));
//	e->a = atof(pick_option(&c, &v, "a", "0"));
//	e->b = atof(pick_option(&c, &v, "b", "0"));
//	e->c = atof(pick_option(&c, &v, "c", "0"));
//	e->R = atof(pick_option(&c, &v, "R", "50"));

	// process input arguments (should be none)
	if (c != 1)
		return fprintf(stderr, "usage:\n\t%s\n", *v);


	// init fonts
	//e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	e->font[0] = reformat_font(*xfont_8x13, UNPACKED);

	// open the window
	struct FTR f = ftr_new_window(2*N, N);
	f.userdata = e;
	f.changed = 1;

	// set event handlers
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "key", event_key);

	return ftr_loop_run(&f) - 1;
}

int main(int c, char *v[]) { return main_vac(c, v); }

// vim:set foldmethod=marker:
