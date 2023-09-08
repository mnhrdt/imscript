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
	// point data (input)
	int n;       // number of points
	float *x;    // point coordinates (x[2*i+0], x[2*i+1]) for 0<=i<n
	//float *c;    // rgb coloring of each point
	int *C;      // label id of each point
	char **f;    // filename/id of each point

	//// computed point data (intermediary)
	//float *y;    // coordinates of inverted points
	//float *z;    // coordinates of the points of the convex hull
	//int m;       // number of points in the convex hull

	//// katz
	//float c[2];  // center of view
	//float r;     // radius of inversion

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
	bool show_debug; // show inverted points and full qhulls

	// ui
	struct bitmap_font font[5]; // from small to large
	int p;                      // whether to show point p
	int px, py;                 // position for showing it
};


// function to reset and center the viewer
static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

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
	e->show_debug = false;

	// x
	e->p = -1;

	f->changed = 1;
}


// xterm-like 16 color palette
static uint8_t palette[16][3] = {
	{0, 0, 0},        //  0  black
	{128, 0, 0},      //  1  dark red
	{0, 128, 0},      //  2  dark green
	{128, 128, 0},    //  3  dark yellow
	{0, 0, 128},      //  4  dark blue
	{128, 0, 128},    //  5  dark magenta
	{0, 128, 128},    //  6  dark cyan
	{192, 192, 192},  //  7  dark white
	{128, 128, 128},  //  8  gray
	{255, 0, 0},      //  9  red
	{0, 255, 0},      // 10  green
	{255, 255, 0},    // 11  yellow
	{0, 0, 255},      // 12  blue
	{255, 0, 255},    // 13  magenta
	{0, 255, 255},    // 14  cyan
	{255, 255, 255},  // 15  white
};


// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 3. algorithms                                                    {{{1

static int compare_points_lexicographically(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	int p = (a[0] > b[0]) - (a[0] < b[0]);
	if (p)
		return p;
	else
		return ((a[1] > b[1]) - (a[1] < b[1]));
}

// oriented area of a triangle
static float det(float a[2], float b[2], float c[2])
{
	float p[2] = {b[0] - a[0], b[1] - a[1]};
	float q[2] = {c[0] - a[0], c[1] - a[1]};
	return p[0] * q[1] - p[1] * q[0];
}

// function to compute the convex hull of a set of points in the plane
// (note: the input points x are sorted in-place)
static int do_the_andrew_parkour(
		float *y,  // output: coordinates of points on the convex hull
		float *x,  // input: list of points in the plane
		int n      // input: number of input points
		)          // return value: number of points in the convex hull
{
	// sort the given points lexicographically (in place)
	qsort(x, n, 2*sizeof*x, compare_points_lexicographically);

	// number of points in the hull found so far
	int r = 0;

	// fill-in the lower hull
	for (int i = 0; i < n; i++)
	{
		while (r >= 2 && det(y+2*(r-2), y+2*(r-1), x+2*i) <= 0)
			r -= 1;
		y[2*r+0] = x[2*i+0];
		y[2*r+1] = x[2*i+1];
		r += 1;
	}

	// fill-in the upper hull
	int k = r + 1; // start of the upper hull
	for (int i = n-2; i >= 0; i--)
	{ //                v--- (only difference with previous loop)
		while (r >= k && det(y+2*(r-2), y+2*(r-1), x+2*i) <= 0)
			r -= 1;
		y[2*r+0] = x[2*i+0];
		y[2*r+1] = x[2*i+1];
		r += 1;
	}

	return r;
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

// draw a segment between two points
static void traverse_circle(int cx, int cy, int r,
		void (*f)(int,int,void*), void *e)
{
	int h = r / sqrt(2);
	for (int i = -h; i <= h; i++)
	{
		int s = sqrt(r*r - i*i);
		f(cx + i, cy + s, e); // upper quadrant
		f(cx + i, cy - s, e); // lower quadrant
		f(cx + s, cy + i, e); // right quadrant
		f(cx - s, cy + i, e); // left quadrant
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

// auxiliary function for drawing a blue pixel
static void plot_pixel_blue(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 255;
	}
}

// auxiliary function for drawing a cyan pixel
static void plot_pixel_cyan(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 255;
		f->rgb[3*idx+2] = 255;
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

// auxiliary function for drawing a gray pixel
static void plot_pixel_gray(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 120;
		f->rgb[3*idx+1] = 120;
		f->rgb[3*idx+2] = 120;
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

// function to draw a blue segment
static void plot_segment_blue(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_blue, f);
}

// function to draw a cyan segment
static void plot_segment_cyan(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_cyan, f);
}

// function to draw a gray segment
static void plot_segment_gray(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_gray, f);
}

// function to draw a green segment
static void plot_segment_green(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_green, f);
}

static void plot_circle_green(struct FTR *f,
		float x, float y, float r)
{
	traverse_circle(x, y, r, plot_pixel_green, f);
}



// Subsection 7.2. Drawing user-interface elements                          {{{2

static void splat_disk(uint8_t *rgb, int w, int h, float p[2], float r,
		uint8_t color[3])
{
	for (int j = -r-1 ; j <= r+1; j++)
	for (int i = -r-1 ; i <= r+1; i++)
	if (hypot(i, j) < r)
	{
		int ii = p[0] + i;
		int jj = p[1] + j;
		if (ii>=0 && jj>=0 && ii<w && jj<h)
		//{
		//	float a = pow(hypot(i, j)/r, 4);
			for (int k = 0; k < 3; k++)
				rgb[3*(w*jj+ii)+k] = color[k];
		//		rgb[3*(w*jj+ii)+k] = a*255 + (1-a)*color[k];
		//}
	}
}


static void draw_red_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	uint8_t red[3] = {255, 0, 0};
	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		map_view_to_window(e, P, e->x + 2*i);
		int Ci = 8 + e->C[i] % 8;
		splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, palette[Ci]);
		//plot_pixel_red(P[0], P[1], f);
		//plot_pixel_red(P[0]+1, P[1], f);
		//plot_pixel_red(P[0]-1, P[1], f);
		//plot_pixel_red(P[0], P[1]+1, f);
		//plot_pixel_red(P[0], P[1]-1, f);
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

static int string_is_filename(char *s)
{
	FILE *f = fopen(s, "r");
	if (!f)
		return 0;
	fclose(f);
	return 1;
}

// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255; // white

	draw_red_points(f);

	if (e->p >= 0 && e->p < e->n)
	{
		char *fn = e->f[e->p];
		if (string_is_filename(fn))
		{
			int w, h;
			uint8_t *a = iio_read_image_uint8(fn, &w, &h);
			for (int j = 0; j < h; j++)
			for (int i = 0; i < w; i++)
			{
				int ii = i + e->px;
				int jj = j + e->py - h;
				if (insideP(f, ii, jj))
				for (int k = 0; k < 3; k++)
					f->rgb[3*(jj*f->w+ii)+k]=a[j*w+i];
			}
			xfree(a);
		} else {
			uint8_t white[3] = {255, 255, 255};
			uint8_t black[3] = {0, 0, 0};
			struct bitmap_font *font = e->font + 4;
			put_string_in_rgb_image(f->rgb, f->w, f->h,
					e->px, e->py - font->height,
					white, black, 0, font, fn);
		}
	}

	//compute_red_points_convex_hull(e);
	//for (int i = 0; i < e->m - 1; i++)
	//{
	//	float P[2], Q[2], Z[4], C[2];
	//	map_view_to_window(e, P, e->z + 2*i);
	//	map_view_to_window(e, Q, e->z + 2*i + 2);
	//	invert_point(Z + 0, e->c, e->r, e->z + 2*i + 0);
	//	invert_point(Z + 2, e->c, e->r, e->z + 2*i + 2);
	//	if (e->show_debug)
	//		plot_segment_gray(f, P[0], P[1], Q[0], Q[1]);
	//	map_view_to_window(e, P, Z + 0);
	//	map_view_to_window(e, Q, Z + 2);
	//	map_view_to_window(e, C, e->c);
	//	if (det(P, Q, C) > 0)
	//	{
	//		plot_segment_blue(f, P[0], P[1], Q[0], Q[1]);
	//		uint8_t blue[3] = {0, 0, 255};
	//		splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, blue);
	//	}
	//	else if (e->show_debug)
	//		plot_segment_cyan(f, P[0], P[1], Q[0], Q[1]);
	//}
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
	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		map_view_to_window(e, P, e->x + 2*i);
		if (hypot(P[0] - x, P[1] - y) < 0.5 + POINT_RADIUS)
			return i;
	}
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
	if (k == 'd') e->show_debug = !e->show_debug;
	//if (k == 'a') e->restrict_to_affine = !e->restrict_to_affine;
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

	int p = hit_point(e, x, y);
	if (p >= 0)
	{
		//fprintf(stderr, "hit point %d (%d %d)\n", p, x, y);
		e->p = p;
		e->px = x;
		e->py = y;
		f->changed = 1;
	}
	if (p < 0 && e->p >= 0)
	{
		e->p = -1;
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

struct data_line {
	float xyrgb[5];
	char name[FILENAME_MAX];
};

static struct data_line *read_data_lines_rgb(FILE *f, int *N)
{
	int n = 0;  // number of lines read
	int m = 8;  // number of lines allocated
	struct data_line *r = xmalloc(m * sizeof*r);
	while (6 == fscanf(f, "%g %g %g %g %g %100[^\n]\n",
				r[n].xyrgb + 0, r[n].xyrgb + 1,
				r[n].xyrgb + 2, r[n].xyrgb + 3, r[n].xyrgb + 4,
				r[n].name)
			)
		if (++n >= m)
			r = xrealloc(r, (m *= 2) * sizeof*r);
	*N = n;
	return r;
}

static struct data_line *read_data_lines_id(FILE *f, int *N)
{
	int n = 0;  // number of lines read
	int m = 8;  // number of lines allocated
	struct data_line *r = xmalloc(m * sizeof*r);
	while (4 == fscanf(f, "%g %g %g %100[^\n]\n",
				r[n].xyrgb + 0, r[n].xyrgb + 1,
				r[n].xyrgb + 2, r[n].name )
			)
		if (++n >= m)
			r = xrealloc(r, (m *= 2) * sizeof*r);
	*N = n;
	return r;
}

int main_cloudette(int argc, char *argv[])
{
	if (argc != 2 && argc != 1) {
		fprintf(stderr, "usage:\n\t%s < xyrgbn.txt\n", *argv);
		return 1;
	}

	// read input data
	int n;
	//struct data_line *l = read_data_lines_rgb(stdin, &n);
	struct data_line *l = read_data_lines_id(stdin, &n);
	fprintf(stderr, "got %d points from stdin\n", n);
	if (n > 10)
	for (int i = 0; i < 10; i++)
		fprintf(stderr, "[%d] %g %g %g \"%s\"\n",
				i,
				l[i].xyrgb[0],
				l[i].xyrgb[1],
				l[i].xyrgb[2],
				l[i].name);

	// initialize state with the given points
	struct viewer_state e[1];
	e->n = n;
	e->x = xmalloc(2*n * sizeof*e->x);
	//e->c = xmalloc(3*n * sizeof*e->c);
	e->C = xmalloc(1*n * sizeof*e->C);
	e->f = xmalloc(1*n * sizeof*e->f);
	for (int i = 0; i < n; i++)
	{
		e->x[2*i+0] = l[i].xyrgb[0];
		e->x[2*i+1] = l[i].xyrgb[1];
		//e->c[3*i+0] = l[i].xyrgb[2];
		//e->c[3*i+1] = l[i].xyrgb[3];
		//e->c[3*i+2] = l[i].xyrgb[4];
		e->C[i] = l[i].xyrgb[2];
		e->f[i] = l[i].name;
	}

	// init fonts
	e->font[0] = reformat_font(*xfont_4x6, UNPACKED);
	e->font[1] = reformat_font(*xfont_6x12, UNPACKED);
	e->font[2] = reformat_font(*xfont_7x13, UNPACKED);
	e->font[3] = reformat_font(*xfont_9x15, UNPACKED);
	e->font[4] = reformat_font(*xfont_10x20, UNPACKED);

	// open the window
	struct FTR f = ftr_new_window(800,600);
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

int main(int c, char *v[]) { return main_cloudette(c, v); }

// vim:set foldmethod=marker:
