// fauxfilet: a simulator of parabolic folds and their dips
// compilation: run "make bin/fauxfilet" in the imscript root


// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

// user interface library
#include "ftr.h"

// image input/output
#include "iio.h"

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

	// cylinder radius
	float R;

	// ui
	struct bitmap_font font[1];
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
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}




// SECTION 3. algorithms                                                    {{{1

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
static float stratum(float x)
{
	return stratum_amfm(x);
}


// SECTION 4. Coordinate Conversions                                        {{{1

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

// take a point of the unfolded cylinder (rectangle), and map it to 3D
static void map_cyl_unfold(struct viewer_state *e, float xyz[3], float ij[2])
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


// SECTION 7. Drawing                                                       {{{1


// Subsection 7.2. Drawing user-interface elements                          {{{2


// classical ``dirt'' palette, as used by dip pickers
static void get_dirt(uint8_t rgb[3], int g)
{
	rgb[1] = 255 - g;
	rgb[0] = g > 127 ? 510 - 2*g : 255;
	rgb[2] = g > 127 ? 0 : 255 - 2*g;
}

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

	// window 0: unfolded cylindrical dip
	// 360x720 starting at 0,0
	for (int j = 0; j < 720; j++)
	for (int i = 0; i < 360; i++)
	{
		float ij[2] = {i, j-360};
		float xyz[3];
		map_cyl_unfold(e, xyz, ij);
		float h = xyz[2] + e->p * xyz[0] * xyz[0];
		float c = stratum(e->f * h);
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
		float h = z + e->p * x*x;
		float c = stratum(e->f * h);
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
	snprintf(buf, 0x200, "f=%g\np=%g\na=%g\nb=%g\nc=%g\nR=%g",
			e->f, e->p, e->a, e->b, e->c, e->R);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			360+0, 0+0, fg, bg, 0, e->font, buf);
}



// SECTION 8. User-Interface Actions and Events                             {{{1


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

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_cloudette_%d.png", c);
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
	fprintf(stderr, "key k=%d m=%d\n", k, m);

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

	// f, p, a, b, c, R hitboxes of font height
	// 0  1  2  3  4  5
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) scale_strata_frequency(e, 1/1.3);
		if (Y == 1) scale_fold_parameter(e, 1/1.3);
		if (Y == 2) shift_angle_a(e, -10);
		if (Y == 3) shift_angle_b(e, -2);
		if (Y == 4) shift_angle_c(e, -2);
		if (Y == 5) scale_radius(e, 1/1.3);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) scale_strata_frequency(e, 1.3);
		if (Y == 1) scale_fold_parameter(e, 1.3);
		if (Y == 2) shift_angle_a(e, 10);
		if (Y == 3) shift_angle_b(e, 2);
		if (Y == 4) shift_angle_c(e, 2);
		if (Y == 5) scale_radius(e, 1.3);
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


// SECTION 9. Image processing


// SECTION 10. Main Program                                                 {{{1

int main_fauxfilet(int argc, char *argv[])
{
	// process input arguments
	if (argc != 1)
		return fprintf(stderr, "usage:\n\t%s\n", *argv);

	// initialize state
	struct viewer_state e[1];
	center_state(e);

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

int main(int c, char *v[]) { return main_fauxfilet(c, v); }

// vim:set foldmethod=marker:
