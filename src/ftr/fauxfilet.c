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
#include "fontu.c" // todo: cherry-pick the required fontu functions
#include "fonts/xfonts_all.c"

// radius of the points
#define POINT_RADIUS 2.3

// zoom factor for zoom-in and zoom-out
#define ZOOM_FACTOR 1.43

// radius scaling factor for the inversion circle
#define RADIUS_FACTOR 1.13

// data structure to store the state of the viewer
struct viewer_state {
	// bullseye parameters
	float f;          // frequency of the strata before folding
	float p;          // parameter of parabolic fold
	float a, b, c;    // euler angles of parabolic sheaf

	// cylinder radius
	float R;

	// window viewport
	float offset[2];  // horizontal (T) and vertical (R) offset
	float scale;      // vertical scale

	// below is unused by now

	// dragging state
	bool dragging_window_point;
	bool dragging_image_point;
	bool dragging_background;
	int dragged_point;
	int drag_handle[2];

	// display options
	int interpolation_order; // 0=nearest, 1=linear, 2=bilinear, 3=bicubic
	bool tile_plane;
	bool show_horizon;
	bool show_grid_points;
	bool restrict_to_affine;
	bool show_debug; // show inverted points and full qhulls

	// ui
	struct bitmap_font font[5]; // from small to large
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

	// viewport
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;

	// drag state
	e->dragging_window_point = false;
	e->dragging_image_point = false;
	e->dragged_point = -1;
	e->dragging_background = false;

	// visualization options
	e->interpolation_order = 0;
	e->tile_plane = false;
	e->show_horizon = false;
	e->show_grid_points = false;
	e->restrict_to_affine = false;
	e->show_debug = false;
}

static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	center_state(e);
	f->changed = 1;
}



// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 3. algorithms                                                    {{{1

static float stratum(float x)
{
	return 127 + 107 * sin(x);
}


// SECTION 4. Coordinate Conversions                                        {{{1

// "view"   : coordinates in the infinite plane where the points are located
// "window" : coordinates in the window, which is a rectangluar piece of "view"
//

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
	//xyz[0] = XYZ[0];
	//xyz[1] = XYZ[1];
	//xyz[2] = XYZ[2];
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
		get_dirt(rgb, g);
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
		if (*rth < e->R)
			get_dirt(rgb, g);
		f->rgb[(j*f->w + 360+i)*3 + 0] = rgb[0];
		f->rgb[(j*f->w + 360+i)*3 + 1] = rgb[1];
		f->rgb[(j*f->w + 360+i)*3 + 2] = rgb[2];
	}

	// hud
	uint8_t fg[3] = {0, 255, 0};
	uint8_t bg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "f=%g p=%g\na=%g b=%g c=%g\nR=%g",
			e->f, e->p, e->a, e->b, e->c, e->R);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			360+0, 0+0, fg, bg, 0, &e->font[4], buf);
}



// SECTION 8. User-Interface Actions and Events                             {{{1


// action: viewport translation
static void change_view_offset(struct viewer_state *e, float dx, float dy)
{
	e->offset[0] += dx;
	e->offset[1] += dy;
}

// action: viewport zoom
static void change_view_scale(struct viewer_state *e, int x, int y, float fac)
{
	float center[2], X[2] = {x, y};
	map_window_to_view(e, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "zoom changed %g\n", e->scale);
}

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


// test whether (x,y) hits some point
static int hit_point(struct viewer_state *e, float x, float y)
{
	//for (int i = 0; i < e->n; i++)
	//{
	//	float P[2];
	//	map_view_to_window(e, P, e->x + 2*i);
	//	if (hypot(P[0] - x, P[1] - y) < 0.5 + POINT_RADIUS)
	//		return i;
	//}
	return -1;
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

	//if (k == 'c') center_view(f);
	if (k == 'f') scale_strata_frequency(e, 1.3);
	if (k == 'F') scale_strata_frequency(e, 1/1.3);
	if (k == 'p') scale_fold_parameter(e, 1.3);
	if (k == 'P') scale_fold_parameter(e, 1/1.3);
	if (k == 'r') scale_radius(e, 1.3);
	if (k == 'R') scale_radius(e, 1/1.3);
	if (k == 'a') shift_angle_a(e, 2);
	if (k == 'A') shift_angle_a(e, -2);
	if (k == 'b') shift_angle_b(e, 2);
	if (k == 'B') shift_angle_b(e, -2);
	if (k == 'c') shift_angle_c(e, 2);
	if (k == 'C') shift_angle_c(e, -2);
	//if (k == 'j') change_view_offset(e, 0, -10);
	//if (k == 'k') change_view_offset(e, 0, 10);
	//if (k == 'h') change_view_offset(e, 10, 0);
	//if (k == 'l') change_view_offset(e, -10, 0);
	//if (k == FTR_KEY_DOWN ) change_view_offset(e, 0, -100);
	//if (k == FTR_KEY_UP   ) change_view_offset(e, 0, 100);
	//if (k == FTR_KEY_RIGHT) change_view_offset(e, -100, 0);
	//if (k == FTR_KEY_LEFT)  change_view_offset(e, 100, 0);
	//if (k == '+') change_view_scale(e, f->w/2, f->h/2, ZOOM_FACTOR);
	//if (k == '-') change_view_scale(e, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);
	//if (k == 'p') e->tile_plane = !e->tile_plane;
	//if (k == 'w') e->show_horizon = !e->show_horizon;
	//if (k >= '0' && k <= '9') e->interpolation_order = k - '0';
	//if (k == '.') e->show_grid_points = !e->show_grid_points;
	//if (k == 'd') e->show_debug = !e->show_debug;
	////if (k == 'a') e->restrict_to_affine = !e->restrict_to_affine;
	if (k == ',') action_screenshot(f);

	e->dragging_window_point = e->dragging_image_point = false;
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

	//// begin dragging a control point in the WINDOW DOMAIN
	//if (k == FTR_BUTTON_LEFT && p >= 0)
	//{
	//	e->dragged_point = p;
	//	e->dragging_window_point = true;
	//}

	//// end dragging a control point in the WINDOW DOMAIN
	//if (e->dragging_window_point && k == -FTR_BUTTON_LEFT)
	//{
	//	drag_point_in_window_domain(e, x, y);
	//	e->dragging_window_point = false;
	//	e->dragged_point = -1;
	//}

	//// begin dragging a control point in the IMAGE DOMAIN
	//if (k == FTR_BUTTON_RIGHT && p >= 0)
	//{
	//	e->dragged_point = p;
	//	e->dragging_image_point = true;
	//}

	// begin dragging a the WINDOW BACKGROUND
	if (k == FTR_BUTTON_LEFT)// && hit_point(e, x, y) < 0)
	{
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		e->dragging_background = true;
	}

	// end dragging the WINDOW BACLGROUND
	if (e->dragging_background && k == -FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->dragging_background = false;
	}

	// radius in/out (if hit), zoom in/out (if no hit)
	if (k == FTR_BUTTON_DOWN)
	{
		change_view_scale(e, x, y, ZOOM_FACTOR);
		//if (hit_point(e, x, y)<0)
		//	change_view_scale(e, x, y, ZOOM_FACTOR);
		//else
		//	change_radius(e, x, y, RADIUS_FACTOR);
	}
	if (k == FTR_BUTTON_UP)
	{
		change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);
		//if (hit_point(e, x, y)<0)
		//	change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);
		//else
		//	change_radius(e, x, y, 1.0/RADIUS_FACTOR);
	}

	f->changed = 1;
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;


	//// drag WINDOW DOMAIN control point (realtime feedback)
	//if (e->dragging_window_point && m & FTR_BUTTON_LEFT)
	//{
	//	drag_point_in_window_domain(e, x, y);
	//	f->changed = 1;
	//}

	// drag WINDOW DOMAIN background (realtime feedback)
	if (e->dragging_background && m & FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		f->changed = 1;
		return;
	}
	f->changed =  1;

//	int p = hit_point(e, x, y);
//	if (p >= 0)
//	{
//		//fprintf(stderr, "hit point %d (%d %d)\n", p, x, y);
//		e->p = p;
//		//e->px = x;
//		//e->py = y;
//		f->changed = 1;
//	}
//	if (p < 0 && e->p >= 0)
//	{
//		e->p = -1;
//		f->changed = 1;
//	}

}

// expose handler
static void event_expose(struct FTR *f, int b, int m, int x, int y)
{
	if (f->changed)
		paint_state(f);
}


// SECTION 9. Image processing

#include "xmalloc.c"
//
//#define PYRAMID_LEVELS 20
//
//struct image_pyramid {
//	float *x[PYRAMID_LEVELS];
//	int w[PYRAMID_LEVELS];
//	int h[PYRAMID_LEVELS];
//	int pd;
//};
//
//static void do_pyramid(struct image_pyramid *p, float *x, int w, int h, int pd)
//{
//	fprintf(stderr, "building a multiscale pyramid %d x %d x %d\n",
//			w, h, pd);
//
//	int nw = 1 + ceil(log2(w+1));
//	int nh = 1 + ceil(log2(h+1));
//	p->n = (nw + nh + abs(nw - nh))/2;
//	p->x = xmalloc(n * sizeof*p->x);
//	p->w = xmalloc(n * sizeof*p->w);
//	p->h = xmalloc(n * sizeof*p->h);
//
//	p->w[0] = w;
//	p->h[0] = h;
//	for (int i = 1; i < p->n; i++)
//
//
//	fprintf(stderr, "\t%d x %d\n", nw, nh);
//}
//
//static void free_pyramid(struct image_pyramid *p)
//{
//	for (int i = 0; i < pyra i++)
//		free(p->x[i]);
//	free(p->x);
//	free(
//}



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
	e->font[0] = reformat_font(*xfont_4x6, UNPACKED);
	e->font[1] = reformat_font(*xfont_6x12, UNPACKED);
	e->font[2] = reformat_font(*xfont_7x13, UNPACKED);
	e->font[3] = reformat_font(*xfont_9x15, UNPACKED);
	e->font[4] = reformat_font(*xfont_10x20, UNPACKED);

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
