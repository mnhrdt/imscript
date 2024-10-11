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

// fourier stuff

// data structure to store the state of the viewer
struct viewer_state {
	// structural parameters
	int w;  // scene width
	int h;  // scene height
	int o;  // scene orientation 0=horizontal, 1=vertical

	// input texture parameters
	int d;          // random seed
	int g;          // gaussian grain of the texture
	int l;          // number of laplacian steps
	int s[2][3];  // two shifts (x,y,sign)
	float p;        // saturation parameter

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
	e->s[0][0] = 100; e->s[0][1] =  0; e->s[0][2] = 1;
	e->s[1][0] =  50; e->s[1][1] = 87; e->s[1][2] = 1;
	e->p = 1;

	// autocorrelation view
	e->r = 20;
	e->z = 1;
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}


// SECTION 2. algorithms                                                    {{{1





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
		int oi = e->w * e->o;
		int oj = e->h * (1 - e->o);
		f->rgb[((j+oj)*f->w + i+oi)*3 + 0] = c[0];
		f->rgb[((j+oj)*f->w + i+oi)*3 + 1] = c[1];
		f->rgb[((j+oj)*f->w + i+oi)*3 + 2] = c[2];
	}


	// free image data
	free(X); free(vX);
	free(x); free(vx);


	// hud
	uint8_t fg[3] = {0, 0, 0};
	uint8_t bg[3] = {127, 127, 127};
	char buf[0x200];
	snprintf(buf, 0x200,
			"d = %d\n"      "g = %g\n"      "l = %d\n"
			"s1 = %d %d\n"  "s2 = %d %d\n"  "p = %g\n"
			"r = %g\n"      "z = %g\n"
			e->d, e->g, e->l, e->

			e->f, e->p, e->a, e->b, e->c, e->R, e->s);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			360+0, 0+0, fg, bg, 0, e->font, buf);
}



// SECTION 5. User-Interface Actions and Events                             {{{1


// action: scale strata frequency
static void scale_strata_frequency(struct viewer_state *e, float f)
{
	e->f *= f;
	//fprintf(stderr, "f = %g\n", e->f);
}
static void scale_fold_parameter(struct viewer_state *e, float f)
{
	e->p *= f;
	//fprintf(stderr, "p = %g\n", e->p);
}
static void shift_angle_a(struct viewer_state *e, float s) { e->a += s; }
static void shift_angle_b(struct viewer_state *e, float s) { e->b += s; }
static void shift_angle_c(struct viewer_state *e, float s) { e->c += s; }
static void scale_radius(struct viewer_state *e, float f) { e->R *= f; }
static void shift_shift_s(struct viewer_state *e, float s) { e->s += s; }
static void shift_shift_z(struct viewer_state *e, float s) { e->z0 += s; }

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
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'z') center_view(f);
	if (k == 'f') scale_strata_frequency(e, 1.3);
	if (k == 'F') scale_strata_frequency(e, 1/1.3);
	if (k == 'p') scale_fold_parameter(e, 1.3);
	if (k == 'P') scale_fold_parameter(e, 1/1.3);
//	if (k == 's') shift_shift_s(e, 1);
//	if (k == 'S') shift_shift_s(e, -1);
	if (k == 'r') scale_radius(e, 1.3);
	if (k == 'R') scale_radius(e, 1/1.3);
	if (k == 'a') shift_angle_a(e, 10);
	if (k == 'A') shift_angle_a(e, -10);
	if (k == 'b') shift_angle_b(e, 2);
	if (k == 'B') shift_angle_b(e, -2);
	if (k == 'c') shift_angle_c(e, 2);
	if (k == 'C') shift_angle_c(e, -2);
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

	// f, p, a, b, c, R, s hitboxes of font height
	// 0  1  2  3  4  5  6
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) scale_strata_frequency(e, 1/1.3);
		if (Y == 1) scale_fold_parameter(e, 1/1.3);
		if (Y == 2) shift_angle_a(e, -10);
		if (Y == 3) shift_angle_b(e, -2);
		if (Y == 4) shift_angle_c(e, -2);
		if (Y == 5) scale_radius(e, 1/1.3);
		if (Y == 6) shift_shift_s(e, -5);
//		if (Y == 7) shift_shift_z(e, -5);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) scale_strata_frequency(e, 1.3);
		if (Y == 1) scale_fold_parameter(e, 1.3);
		if (Y == 2) shift_angle_a(e, 10);
		if (Y == 3) shift_angle_b(e, 2);
		if (Y == 4) shift_angle_c(e, 2);
		if (Y == 5) scale_radius(e, 1.3);
		if (Y == 6) shift_shift_s(e, 5);
//		if (Y == 7) shift_shift_z(e, 5);
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

int main_vac_noninteractive(int c, char *v[])
{
	// initialize state (sets dummy default arguments)
	struct viewer_state e[1];
	center_state(e);
	e->stratum = 0;

	// process named arguments
	int w = atoi(pick_option(&c, &v, "w", "360"));
	int h = atoi(pick_option(&c, &v, "h", "720"));
	char *filename_out = pick_option(&c, &v, "o", "-");
	char *filename_stratum = pick_option(&c, &v, "S", "");
	if (*filename_stratum)
	{
		uint8_t *iio_read_image_uint8(char*,int*,int*);
		e->stratum = iio_read_image_uint8(filename_stratum,
				&e->stratum_w, &e->stratum_h);
		e->f = 1;
	}
	e->f = atof(pick_option(&c, &v, "f", "0.1"));
	e->p = atof(pick_option(&c, &v, "p", "0.001"));
	e->a = atof(pick_option(&c, &v, "a", "0"));
	e->b = atof(pick_option(&c, &v, "b", "0"));
	e->c = atof(pick_option(&c, &v, "c", "0"));
	e->R = atof(pick_option(&c, &v, "R", "50"));

	// fill-in output image
	uint8_t *x = malloc(3*w*h);
	paint_cylinder(x, w, h, e);

	// save output image
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(filename_out, x, w, h, 3);

	return 0;
}

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

	// if -n option, run noninteractively
	if (pick_option(&c, &v, "n", 0))
		return main_vac_noninteractive(c, v);

	// initialize state
	struct viewer_state e[1];
	center_state(e);
	e->stratum = 0;

	// extract named arguments
	char *filename_out = pick_option(&c, &v, "o", "-");
	char *filename_stratum = pick_option(&c, &v, "S", "");
	if (*filename_stratum)
	{
		uint8_t *iio_read_image_uint8(char*,int*,int*);
		e->stratum = iio_read_image_uint8(filename_stratum,
				&e->stratum_w, &e->stratum_h);
		e->f = 1;
	}
	e->f = atof(pick_option(&c, &v, "f", "0.1"));
	e->p = atof(pick_option(&c, &v, "p", "0.001"));
	e->a = atof(pick_option(&c, &v, "a", "0"));
	e->b = atof(pick_option(&c, &v, "b", "0"));
	e->c = atof(pick_option(&c, &v, "c", "0"));
	e->R = atof(pick_option(&c, &v, "R", "50"));

	// process input arguments (should be none)
	if (c != 1)
		return fprintf(stderr, "usage:\n\t%s\n", *v);


	// init fonts
	e->font[0] = reformat_font(*xfont_10x20, UNPACKED);

	// open the window
	struct FTR f = ftr_new_window(1600,800);
	f.userdata = e;
	f.changed = 1;

	// set event handlers
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "key", event_key);

	return ftr_loop_run(&f);
}

int main(int c, char *v[]) { return main_vac(c, v); }

// vim:set foldmethod=marker:
