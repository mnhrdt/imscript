// vpoints: visualize a D-dimensional set of N points


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

#include "iio.h"

// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


#include "random.c"

//#define OMIT_BLUR_MAIN
//#include "blur.c"
//
//#define OMIT_PPSMOOTH_MAIN
//#include "ppsmooth.c"



// data structure to store the state of the viewer
struct viewer_state {
	// structural parameters
	int w;    // scene width
	int h;    // scene height

	// data
	int d;    // dimension of the space
	int n;    // number of points
	float *x; // point coordinates (array nxd)

	// visualizer state
	int   r;  // random seed
	float *c; // center
	float λ;  // scale
	float *p; // first projection direction
	float *q; // second projection direction

	// ui
	int show_basis;
	float point_radius;
	struct bitmap_font font[1];
	float offset[2], scale;
};


// function to reset and center the viewer
static void center_state(struct viewer_state *e)
{
	e->offset[0] = e->w / 2;
	e->offset[1] = e->h / 2;
	e->scale = fmin(e->w, e->h) / 3;

	e->r = 1;

	for (int i = 0; i < e->d; i++)
		e->c[i] = 0;
	e->λ = 1;

	e->show_basis = 1;
	e->point_radius = 2.3;
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}


// SECTION 2. algorithms                                                    {{{1

static float vnorm(float *x, int d)
{
	return d ? hypot(*x, vnorm(x+1,d-1)) : 0;
}

static float sprod(float *x, float *y, int d)
{
	return d?*x**y+sprod(x+1,y+1,d-1):0;
}

static void create_random_directions(struct viewer_state *e)
{
	xsrand(e->r);
	for (int i = 0; i < e->d; i++) e->p[i] = random_normal();
	for (int i = 0; i < e->d; i++) e->q[i] = random_normal();
	float np = vnorm(e->p, e->d);
	for (int i = 0; i < e->d; i++) e->p[i] /= np;
	float pq = sprod(e->p, e->q, e->d);
	for (int i = 0; i < e->d; i++) e->q[i] -= pq * e->p[i];
	float nq = vnorm(e->q, e->d);
	for (int i = 0; i < e->d; i++) e->q[i] /= nq;
	np = vnorm(e->p, e->d);
	nq = vnorm(e->q, e->d);
	pq = sprod(e->p, e->q, e->d);
	assert(fabs(np - 1) < 0.00001);
	assert(fabs(nq - 1) < 0.00001);
	assert(fabs(pq - 0) < 0.00001);
}

// change from plane coordinates to window coordinates
static void map_view_to_window(struct viewer_state *e, float y[2], float x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = e->offset[k] + e->scale * x[k];
}

// change from window coordinates to plane coordinates
static void map_window_to_view(struct viewer_state *e, float y[2], float x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = ( x[k] - e->offset[k] ) / e->scale;
}


// Subsection 7.2. Drawing user-interface elements                          {{{2
//
static void splat_disk(uint8_t *rgb, int w, int h, float p[2], float r,
		uint8_t color[3])
{
	for (int j = -r-1 ; j <= r+1; j++)
	for (int i = -r-1 ; i <= r+1; i++)
	if (hypot(i, j) < r)
	{
		int ii = p[0] + i;
		int jj = p[1] + j;
		float R = hypot(ii - p[0], jj - p[1]);
		//if (R >= r) continue;
		if (ii>=0 && jj>=0 && ii<w && jj<h)
		{
			//float a = pow(R/r, 2);
			for (int k = 0; k < 3; k++)
				rgb[3*(w*jj+ii)+k] = color[k];
		//		rgb[3*(w*jj+ii)+k] = a*255 + (1-a)*color[k];
		//		rgb[3*(w*jj+ii)+k] = a*rgb[3*(w*jj+ii)+k] + (1-a)*color[k];
		}
	}
}



// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas (white)
	for (int i = 0 ; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 255;

	create_random_directions(e);

	// draw unit sphere in light gray
	uint8_t gray[3] = {200, 200, 200};
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		float P[2], ij[2] = {i, j};
		map_window_to_view(e, P, ij);
		if (hypot(P[0], P[1]) < 1)
		{
			f->rgb[3*(j*f->w+i)+0] = gray[0];
			f->rgb[3*(j*f->w+i)+1] = gray[1];
			f->rgb[3*(j*f->w+i)+2] = gray[2];
		}
	}


	// draw points in red
	uint8_t red[3] = {255, 0, 0};
	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		P[0] = sprod(e->p, e->x + i*e->d, e->d);
		P[1] = sprod(e->q, e->x + i*e->d, e->d);

		float ij[2];
		map_view_to_window(e, ij, P);
		splat_disk(f->rgb, f->w, f->h, ij, e->point_radius, red);
	}

	// draw the canonical basis in black
	uint8_t black[3] = {0, 0, 0};
	for (int i = 0; i < e->d; i++)
	{
		float P[2] = { e->p[i], e->q[i]} ;

		float ij[2];
		map_view_to_window(e, ij, P);
		splat_disk(f->rgb, f->w, f->h, ij, 2*e->point_radius, black);
	}


	//// space for images and their visualizations
	//float *x = xmalloc(e->w * e->h * sizeof*x);
	//uint8_t *vx = xmalloc(3* e->w * e->h);

	//// perform all the computaitons
	//fill_base_texture(x, e);
	////distort_base_texture(x, e);
	//display_base_texture(vx, x, e);

	//// window 0: base texture
	//// WxH starting at (0,0)
	//for (int j = 0; j < e->h; j++)
	//for (int i = 0; i < e->w; i++)
	//{
	//	uint8_t *c = vx + 3*(j*e->w + i);
	//	f->rgb[(j*f->w + i)*3 + 0] = c[0];
	//	f->rgb[(j*f->w + i)*3 + 1] = c[1];
	//	f->rgb[(j*f->w + i)*3 + 2] = c[2];
	//}


	//// free image data
	//free(x); free(vx);


	//// hud
	////uint8_t fg[3] = {0, 0, 0};
	////uint8_t bg[3] = {127, 127, 127};
	//uint8_t fg[3] = {0, 255, 0};
	////uint8_t bg[3] = {100, 100, 100};
	//uint8_t bg[3] = {0, 0, 0};
	////uint8_t fg[3] = {255, 0, 255};
	////uint8_t *bg = 0;
	//char buf[0x200];
	//snprintf(buf, 0x200,
	//		"d (random seed) = %d\n"
	//		"a (stability)   = %g\n"
	//		"b (skewness)    = %g\n"
	//		"k (kernel)      = %s\n"
	//		"s (scale)       = %g\n"
	//		"p (palette)     = %d\n"
	//		"r (saturation)  = %g\n",
	//		e->d, e->a, e->b,
	//		global_kernel_names[e->k],
	//		e->s, e->p, e->r);
	//put_string_in_rgb_image(f->rgb, f->w, f->h,
	//		0+0, 0+0, fg, bg, 0, e->font, buf);
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

static void shift_random_seed(struct viewer_state *e, float s) { e->r += s; }
//static void shift_param_b(struct viewer_state *e, float s)
//{
//	e->b += s;
//	e->b = fmin(e->b, 1);
//	e->b = fmax(e->b, -1);
//}
//static void shift_param_a(struct viewer_state *e, float s)
//{
//	float C = 20;
//	float x = -(C-2) + (2*C-2)/e->a;
//	x += s;
//	x = fmax(x, 1);
//	e->a = (2*C-2)/(C-2 + x);
//}
//static void shift_param_k(struct viewer_state *e, float s)
//{
//	e->k += s;
//	if (e->k < 0) e->k = 0;
//	if (e->k >= NUM_KERNELS) e->k = NUM_KERNELS - 1;
//}
//static void toggle_palette(struct viewer_state *e, int s)
//{
//	e->p += s;
//}
//static void shift_saturation(struct viewer_state *e, float s)
//{
//	e->r += s;
//	if (e->k < 0) e->k = 0;
//	if (e->k > 49) e->k = 49;
//}
//static void shift_num_laplacians(struct viewer_state *e, int s) { e->l += s; }
//static void shift_mask_radius(struct viewer_state *e, int s) { e->r += s; }
//static void shift_quantization_q(struct viewer_state *e, int s) { e->q += s; }
//static void scale_gaussian_grain(struct viewer_state *e, float f) { e->g *= f; }
//static void scale_param_s(struct viewer_state *e, float f) { e->s *= f; }
//static void scale_saturation_p(struct viewer_state *e, float f) { e->p *= f; }
//static void scale_saturation_z(struct viewer_state *e, float f) { e->z *= f; }
//static void scale_fractional_a(struct viewer_state *e, float f) { e->a *= f; }

//static void shift_shift(struct viewer_state *e, int p, int q, int s)
//{
//	if (p < 0) p = 0;
//	if (p > 1) p = 1;
//	if (q < 0) q = 0;
//	if (q > 2) q = 2;
//	int f = q < 2 ? 4 : 1;
//	e->s[p][q] += f*s;
//}

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	int p = getpid();
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_vstab_%d_%d.png", p, c);
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
	int X = x / e->font->width;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) shift_random_seed(e, -1);
//		if (Y == 1) shift_param_a(e, -1);
//		if (Y == 2) shift_param_b(e, -0.1);
//		if (Y == 3) shift_param_k(e, -1);
//		if (Y == 4) scale_param_s(e, 1/1.3);
//		if (Y == 5) toggle_palette(e, -1);
//		if (Y == 6) shift_saturation(e, -0.125);
//		if (Y == 1) scale_gaussian_grain(e, 1/1.3);
//		if (Y == 2) shift_num_laplacians(e, -1);
//		if (Y == 4) shift_shift(e, 1, (X-20)/4, -1);
//		if (Y == 5) scale_saturation_p(e, 1/1.3);
//		if (Y == 6) shift_quantization_q(e, -1);
//		if (Y == 7) scale_fractional_a(e, 1/1.1);
//		if (Y == 8) shift_mask_radius(e, -1);
//		if (Y == 9) scale_saturation_z(e, 1/1.3);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) shift_random_seed(e, 1);
//		if (Y == 1) shift_param_a(e, 1);
//		if (Y == 2) shift_param_b(e, 0.1);
//		if (Y == 3) shift_param_k(e, 1);
//		if (Y == 4) scale_param_s(e, 1.3);
//		if (Y == 5) toggle_palette(e, 1);
//		if (Y == 6) shift_saturation(e, 0.125);
//		if (Y == 1) scale_gaussian_grain(e, 1.3);
//		if (Y == 2) shift_num_laplacians(e, 1);
//		if (Y == 3) shift_shift(e, 0, (X-20)/4, 1);
//		if (Y == 4) shift_shift(e, 1, (X-20)/4, 1);
//		if (Y == 5) scale_saturation_p(e, 1.3);
//		if (Y == 6) shift_quantization_q(e, 1);
//		if (Y == 7) scale_fractional_a(e, 1.1);
//		if (Y == 8) shift_mask_radius(e, 1);
//		if (Y == 9) scale_saturation_z(e, 1.3);
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

static char *help_string_name     = "vpoints";
static char *help_string_version  = "vpoints 1.0\n\nWritten by mnhrdt";
static char *help_string_oneliner = "visualize a n-dimensional point cloud";
static char *help_string_usage    = "usage:\n\tvpoints < points";
static char *help_string_long     =
"Vstab is an interface for blurred and quantified stable noise.\n"
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
int main_vpoints(int c, char *v[])
{
	// if requested, print help
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

//	// if -n option, run noninteractively
//	if (pick_option(&c, &v, "n", 0))
//		return main_vac_noninteractive(c, v);

	// initialize state
	struct viewer_state e[1];
	e->x = iio_read_image_float("-", &e->d, &e->n);
	e->c = xmalloc(e->d * sizeof*e->c);
	e->p = xmalloc(e->d * sizeof*e->p);
	e->q = xmalloc(e->d * sizeof*e->q);

	int W = 800;
	e->w = W;
	e->h = W;
	center_state(e);

	// extract named arguments
	//char *filename_out = pick_option(&c, &v, "o", "-");
	//char *filename_base = pick_option(&c, &v, "S", "");
	//if (*filename_base)
	//	e->stratum = iio_read_image_float(filename_base, &e->w, &e->h);
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
	//e->font[0] = reformat_font(*xfont_7x13, UNPACKED);

	// open the window
	struct FTR f = ftr_new_window(e->w, e->h);
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

int main(int c, char *v[]) { return main_vpoints(c, v); }

// vim:set foldmethod=marker:
