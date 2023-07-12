// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//// user interface library
#include "ftr.h"

// radius of the disks that are displayed around control points
#define DISK_RADIUS 7

// zoom factor for zoom-in and zoom-out
#define ZOOM_FACTOR 1.43


// data structure to store the state of the viewer
struct viewer_state {
	// point data
	int n;       // number of points
	float *x;    // point coordinates (x[2*i+0], x[2*i+1]) for 0<=i<n

	// katz
	float c[2];  // center of view
	float r;     // radius of inversion

	// window viewport
	float offset[2];
	float scale;

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
};


// function to reset and center the viewer
static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// katz
	e->c[0] = 100;
	e->c[1] = 100;
	e->r = 100;

	// drag state
	e->dragging_window_point = false;
	e->dragging_image_point = false;
	e->dragged_point = -1;
	e->dragging_background = false;

	// viewport
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;

	// visualization options
	e->interpolation_order = 0;
	e->tile_plane = false;
	e->show_horizon = false;
	e->show_grid_points = false;
	e->restrict_to_affine = false;

	f->changed = 1;
}


// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 2. Linear algebra                                                {{{1


// compute the vector product of two vectors
static void vector_product(float axb[3], float a[3], float b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}



// SECTION 3. Coordinate Conversions                                        {{{1

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


// SECTION 7. Drawing                                                       {{{1

// Subsection 7.1. Drawing segments                                         {{{2

// generic function to traverse a segment between two pixels
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (qx - px > qy - py || px - qx > qy - py) { // horizontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px + j*slope), j+py, e);
		}
	}
}

// auxiliary function for drawing a red pixel
static void plot_pixel_red(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}

// auxiliary function for drawing a green pixel
static void plot_pixel_green(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 128;
		f->rgb[3*idx+2] = 0;
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

// function to draw a green segment
static void plot_segment_green(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_green, f);
}



// Subsection 7.2. Drawing user-interface elements                          {{{2

// draw the view point
static void draw_view_center(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	float P[2];
	map_view_to_window(e, P, e->c);

	// grey circle
	int side = DISK_RADIUS;
	for (int j = -side-1 ; j <= side+1; j++)
	for (int i = -side-1 ; i <= side+1; i++)
	if (hypot(i, j) < side)
	{
		int ii = P[0] + i;
		int jj = P[1] + j;
		if (insideP(f, ii, jj))
			for (int c = 0; c < 3; c++)
				f->rgb[3*(f->w*jj+ii)+c] = 127;
	}

	// central green dot
	int ii = P[0];
	int jj = P[1];
	if (insideP(f, ii, jj))
		f->rgb[3*(f->w*jj+ii)+1]=255;
}

static void draw_red_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		map_view_to_window(e, P, e->x + 2*i);
		plot_pixel_red(P[0], P[1], f);
		plot_pixel_red(P[0]+1, P[1], f);
		plot_pixel_red(P[0]-1, P[1], f);
		plot_pixel_red(P[0], P[1]+1, f);
		plot_pixel_red(P[0], P[1]-1, f);
	}
}

//// draw four red segments connecting the control points
//static void draw_four_red_segments(struct FTR *f)
//{
//	struct viewer_state *e = f->userdata;
//
//	int o[5] = {0, 1, 3, 2, 0}; // order of the points
//	for (int p = 0; p < 4; p++)
//	{
//		float x0[2], xf[2];
//		map_view_to_window(e, x0, e->c[o[p+0]]);
//		map_view_to_window(e, xf, e->c[o[p+1]]);
//		plot_segment_red(f, x0[0], x0[1], xf[0], xf[1]);
//	}
//}



// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255; // white

	draw_red_points(f);
	draw_view_center(f);
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

// action: drag a point in the window domain
static void drag_point_in_window_domain(struct viewer_state *e, int x, int y)
{
	float X[2] = {x, y};
	map_window_to_view(e, e->c, X);
}

// test whether (x,y) is inside one of the four control disks
static int hit_point(struct viewer_state *e, float x, float y)
{
	float P[2];
	map_view_to_window(e, P, e->c);
	if (hypot(P[0] - x, P[1] - y) < 2 + DISK_RADIUS)
		return 0;
	return -1;
}

// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'c') center_view(f);
	if (k == 'j') change_view_offset(e, 0, -10);
	if (k == 'k') change_view_offset(e, 0, 10);
	if (k == 'h') change_view_offset(e, 10, 0);
	if (k == 'l') change_view_offset(e, -10, 0);
	if (k == FTR_KEY_DOWN ) change_view_offset(e, 0, -100);
	if (k == FTR_KEY_UP   ) change_view_offset(e, 0, 100);
	if (k == FTR_KEY_RIGHT) change_view_offset(e, -100, 0);
	if (k == FTR_KEY_LEFT)  change_view_offset(e, 100, 0);
	if (k == '+') change_view_scale(e, f->w/2, f->h/2, ZOOM_FACTOR);
	if (k == '-') change_view_scale(e, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);
	if (k == 'p') e->tile_plane = !e->tile_plane;
	if (k == 'w') e->show_horizon = !e->show_horizon;
	if (k >= '0' && k <= '9') e->interpolation_order = k - '0';
	if (k == '.') e->show_grid_points = !e->show_grid_points;
	//if (k == 'a') e->restrict_to_affine = !e->restrict_to_affine;

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

	int p = hit_point(e, x, y);

	// begin dragging a control point in the WINDOW DOMAIN
	if (k == FTR_BUTTON_LEFT && p >= 0)
	{
		e->dragged_point = p;
		e->dragging_window_point = true;
	}

	// end dragging a control point in the WINDOW DOMAIN
	if (e->dragging_window_point && k == -FTR_BUTTON_LEFT)
	{
		drag_point_in_window_domain(e, x, y);
		e->dragging_window_point = false;
		e->dragged_point = -1;
	}

	// begin dragging a control point in the IMAGE DOMAIN
	if (k == FTR_BUTTON_RIGHT && p >= 0)
	{
		e->dragged_point = p;
		e->dragging_image_point = true;
	}

	// begin dragging a the WINDOW BACKGROUND
	if (k == FTR_BUTTON_LEFT && hit_point(e, x, y) < 0)
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

	// zoom in/out
	if (k == FTR_BUTTON_DOWN) change_view_scale(e, x, y, ZOOM_FACTOR);
	if (k == FTR_BUTTON_UP) change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);

	f->changed = 1;
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// drag WINDOW DOMAIN control point (realtime feedback)
	if (e->dragging_window_point && m & FTR_BUTTON_LEFT)
	{
		drag_point_in_window_domain(e, x, y);
		f->changed = 1;
	}

	// drag WINDOW DOMAIN background (realtime feedback)
	if (e->dragging_background && m & FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		f->changed = 1;
	}
}

// expose handler
static void event_expose(struct FTR *f, int b, int m, int x, int y)
{
	if (f->changed)
		paint_state(f);
}


// SECTION 9. Image processing

//#include "xmalloc.c"
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

#include "parsenumbers.c"  // read_ascii_floats

// main function
int main_pkatz(int argc, char *argv[])
{
	if (argc != 2 && argc != 1) {
		fprintf(stderr, "usage:\n\t%s < points.txt\n", *argv);
		return 1;
	}

	// initialize state with the given points
	struct viewer_state e[1];
	e->x = read_ascii_floats(stdin, &e->n);
	e->n /= 2;
	fprintf(stderr, "read %d points from stdin\n", e->n);

	// open the window
	struct FTR f = ftr_new_window(512,512);
	f.userdata = e;
	center_view(&f);

	// set event handlers
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "key", event_key);

	return ftr_loop_run(&f);
}

int main(int c, char *v[]) { return main_pkatz(c, v); }

// vim:set foldmethod=marker:
