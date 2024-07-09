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

// data structure to store the state of the viewer
struct viewer_state {
	// bullseye parameters
	float f;          // frequency of the strata before folding
	float p;          // parameter of parabolic fold
	float a, b, c;    // euler angles of parabolic sheaf
	float s;          // horizontal shift of the fold
	float z0;          // vertical shift of the fold

	// cylinder radius
	float R;

	// ui
	struct bitmap_font font[1];

	// data
	int stratum_w;
	int stratum_h;
	uint8_t *stratum;
};


// function to reset and center the viewer
static void center_state(struct viewer_state *e)
{
	// bullseye
	e->f = 0.1;  // width of the strata = 10 pixels
	e->p = 0.0001;    // flat strata (zero parabolic fold)
	e->a = 0;    // euler angles zeroed (vertical cylinter)
	e->b = 0;
	e->c = 0;
	e->R = 20;

	e->s = 0;
	e->z0 = 0;
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}


// SECTION 2. algorithms                                                    {{{1

// parallel parabolic fold of parameter p, horizontal shift s and vert shift z0
static float fold(float x, float y, float z, float p, float s, float z0)
{
	x -= s;
	z -= z0;
	return z + p*x*x;
}

// basic sinusoidal wave
static float stratum_sin(float x)
{
	return 127 + 107 * sin(x);
}

// frequency modulated sinusoid, non-periodic
static float stratum_fm(float x)
{
	return 127 + 107 * sin(x+cos(M_PI*x));
}

// amplitude modulated sinusoid, non-periodic
static float stratum_am(float x)
{
	return 127 + 107 * sin(x)*(0.5+0.5*cos(M_PI*x));
}

// frequency and amplitude modulated
static float stratum_amfm(float x)
{
	return 127 + 107 * sin(1+x+cos(M_PI*x))*(0.5+0.5*cos(M_2_PI*x+sin(x)));
}

// TODO: interactive access to this choice, and to the modulation parameters
static float stratum(struct viewer_state *e, float x)
{
	x *= e->f;
	if (e->stratum)
	{
		float fj = x + e->stratum_h/2;
		int j = fj;
		int i = e->stratum_w/2;
		if (j < 0) return 0; //j = 0;
		if (j >= e->stratum_h-1) return 0;//j = e->stratum_h - 1;
		float a = fj - j;
		//if (a < 1 
		return (1-a)*e->stratum[j*e->stratum_w+i]
			+a*e->stratum[(j+1)*e->stratum_w+i];
	}
	else
		return stratum_amfm(x);
}




// SECTION 3. Coordinate Conversions                                        {{{1

// y = Ax
static void matvec33(float y[3], float A[3][3], float x[3])
{
	for (int i = 0; i < 3; i++)
		y[i] = 0;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		y[i] += A[i][j] * x[j];
}

// apply euler rotations to a 3d point
static void map_euler(struct viewer_state *e, float X[3], float x[3], int dir)
{
	static float ca, sa, cb, sb, cc, sc;

	float a = dir * e->a * M_PI / 180;
	float b = dir * e->b * M_PI / 180;
	float c = dir * e->c * M_PI / 180;

	float R[3][3][3] = {
		{
			{ cos(a), sin(a), 0 },
			{ -sin(a), cos(a), 0 },
			{ 0, 0, 1}
		},
		{
			{ cos(b), 0, -sin(b) },
			{ 0, 1, 0},
			{ sin(b), 0, cos(b) }
		},
		{
			{ 1, 0, 0},
			{ 0, cos(c), sin(c) },
			{ 0, -sin(c), cos(c) }
		}
	};

	float t[3];
	if (dir > 0) {
		matvec33(X, R[0], x);
		matvec33(t, R[1], X);
		matvec33(X, R[2], t);
	} else {
		matvec33(X, R[2], x);
		matvec33(t, R[1], X);
		matvec33(X, R[0], t);
	}
}

// take a point of the unwrapped cylinder (rectangle), and map it to 3D
static void map_cyl_unwrap(struct viewer_state *e, float xyz[3], float ij[2])
{
	float XYZ[3] = {
		e->R * cos(ij[0] * M_PI / 180),
		e->R * sin(ij[0] * M_PI / 180),
		ij[1]
	};
	map_euler(e, xyz, XYZ, 1);
}

// get the cylindrical coordinates of a 3d point
static void map_getcyl(struct viewer_state *e, float rth[3], float xyz[3])
{
	float XYZ[3];
	map_euler(e, XYZ, xyz, -1);
	rth[0] = hypot(XYZ[0], XYZ[1]);
	rth[1] = atan2(XYZ[1], XYZ[0]) * 180 / M_PI;
	rth[2] = XYZ[2];
}


// SECTION 4. Drawing                                                       {{{1

// classical ``dirt'' palette, as used by dip pickers
static void get_dirt(uint8_t rgb[3], int g)
{
	rgb[1] = 255 - g;
	rgb[0] = g > 127 ? 510 - 2*g : 255;
	rgb[2] = g > 127 ? 0 : 255 - 2*g;
}

// paint a cylinder on a rgb image
static void paint_cylinder(uint8_t *x, int w, int h,
		struct viewer_state *e)
{
	// window 0: unwrapped cylindrical dip
	// WxH=360x720 starting at 0,0
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float ij[2] = {i*360.0/w, j-h/2.0};
		float xyz[3];
		map_cyl_unwrap(e, xyz, ij);
		float h = fold(xyz[0], xyz[1], xyz[2], e->p, e->s, e->z0);
		float c = stratum(e, h);
		uint8_t g = c;
		uint8_t rgb[3];
		get_dirt(rgb, 255-g);
		x[(j*w + i)*3 + 0] = rgb[0];
		x[(j*w + i)*3 + 1] = rgb[1];
		x[(j*w + i)*3 + 2] = rgb[2];
	}
}


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

	// window 0: unwrapped cylindrical dip
	// 360x720 starting at 0,0
	for (int j = 0; j < 720; j++)
	for (int i = 0; i < 360; i++)
	{
		float ij[2] = {i, j-360};
		float xyz[3];
		map_cyl_unwrap(e, xyz, ij);
		float h = fold(xyz[0], xyz[1], xyz[2], e->p, e->s, e->z0);
		float c = stratum(e, h);
		uint8_t g = c;
		uint8_t rgb[3];
		get_dirt(rgb, 255-g);
		f->rgb[(j*f->w + i)*3 + 0] = rgb[0];
		f->rgb[(j*f->w + i)*3 + 1] = rgb[1];
		f->rgb[(j*f->w + i)*3 + 2] = rgb[2];
	}

	// window 1: strata
	// 720x720 starting at 360,0
	for (int j = 0; j < 720; j++)
	for (int i = 0; i < 720; i++)
	{
		float x = i - 360;
		float z = j - 360;
		float h = fold(x, 0, z, e->p, e->s, e->z0);
		float c = stratum(e, h);
		uint8_t g = c;
		float xyz[3] = { x, 0, z};
		float rth[3];
		map_getcyl(e, rth, xyz);
		uint8_t rgb[3] = {g, g, g};
		if (*rth < e->R && fabs(rth[2]) < 360)
			get_dirt(rgb, 255-g);
		else {rgb[0]*=0.8;rgb[1]*=0.8;rgb[2]*=0.8;}
		f->rgb[(j*f->w + 360+i)*3 + 0] = rgb[0];
		f->rgb[(j*f->w + 360+i)*3 + 1] = rgb[1];
		f->rgb[(j*f->w + 360+i)*3 + 2] = rgb[2];
	}

	// hud
	uint8_t fg[3] = {0, 255, 0};
	uint8_t bg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "f=%g\np=%g\na=%g\nb=%g\nc=%g\nR=%g\ns=%g\n",
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
static char *help_string_oneliner = "interactive display parabolic folds";
static char *help_string_usage    = "usage:\n\tvac";
static char *help_string_long     =
"Fauxfilet is an interface for exploring parabolic folds.\n"
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
" -w X\tset image width (default=360)\n"
" -h X\tset image height (default=720)\n"
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
	struct FTR f = ftr_new_window(1080,720);
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
