// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

//// user interface library
#include "ftr.h"


// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


// radius of the disks that are displayed around control points
#define DISK_RADIUS 7.3

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

	// general data used for all methods
	float c[2];  // center of view

	int mode; // 0=katz, 1=biasutti

	// katz-specific data
	float r;     // radius of inversion
	float *y;    // coordinates of inverted points
	float *z;    // coordinates of the points of the convex hull
	int m;       // number of points in the convex hull

	// biasutti-specific data
	int N;   // number of knn points
	float alpha; // a-parameter (between 0 and 1)
	float *P; // coordinates of chosen points
	//int *closest; // index of closest point
	//int *farthest; // index of farthest point

	// debug stuff
	int i; // current debug point
	float Ni[100][2]; // coordinates of debug point angle neighbors




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

	struct bitmap_font font[1];
	int negative; // for dark-background mode
};


// function to reset and center the viewer
static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// katz
	e->c[0] = 100;
	e->c[1] = 100;
	e->r = 400;

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

// SECTION 3. Katz algorithm                                                {{{1


static void invert_point_classic(float y[2], float c[2], float R, float x[2])
{
	float r = hypot(x[0] - c[0], x[1] - c[1]);
	y[0] = c[0] + R*R*(x[0] - c[0])/r/r;
	y[1] = c[1] + R*R*(x[1] - c[1])/r/r;
}

static void invert_point_gamma(float y[2], float c[2], float R, float x[2])
{
	// warning: don't use
	float gamma = 0.01;
	float r = hypot(x[0] - c[0], x[1] - c[1]);
	y[0] = c[0] + R*(x[0] - c[0])/pow(r, gamma);
	y[1] = c[1] + R*(x[1] - c[1])/pow(r, gamma);
}

// "katz inversion"
// Katz-Leifman-Tal
// Mesh segmentation using feature point and core extraction
// The Visual Computer (2005)
static void invert_point_flip(float y[2], float c[2], float R, float x[2])
{
	float r = hypot(x[0] - c[0], x[1] - c[1]);
	y[0] = x[0] + 2*(R-r)*(x[0] - c[0])/r;
	y[1] = x[1] + 2*(R-r)*(x[1] - c[1])/r;
}

static void invert_point(float y[2], float c[2], float R, float x[2])
{
	invert_point_flip(y, c, R, x);
	//invert_point_gamma(y, c, R, x);
}


static void invert_points(float *y, float c[2], float R, float *x, int n)
{
	for (int i = 0; i < n; i++)
		invert_point(y + 2*i, c, R, x + 2*i);
}

static void compute_points_inversion(struct viewer_state *e)
{
	invert_points(e->y, e->c, e->r, e->x, e->n);
}


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

static void compute_red_points_convex_hull(struct viewer_state *e)
{
	e->m = do_the_andrew_parkour(e->z, e->y, e->n);
	if (e->show_debug)
	fprintf(stderr, "convex hull of %d points has %d points (%g%%)\n",
			e->n, e->m, e->m * 100.0 / e->n);
}


// SECTION 4. Biasutti's algorithm                                          {{{1

// angle comparison function
static float angle_dist(float a, float b)
{
	return M_PI - fabs(fmod(fabs(a - b), 2*M_PI) - M_PI);
}

static float point_dist(float *a, float *b)
{
	return hypot(a[0] - b[0], a[1] - b[1]);
}

// main algorithm
static void biasutti(struct viewer_state *e)
{
	// compute and sort the angles of each point
	float a[e->n][2]; // [angle, original index]
	for (int i = 0; i < e->n; i++)
	{
		a[i][0] = atan2(e->x[2*i+1] - e->c[1], e->x[2*i+0] - e->c[0]);
		a[i][1] = i;
	}
	assert(sizeof*a == 2*sizeof(float));
	qsort(a, e->n, sizeof*a, compare_points_lexicographically);
	//for (int i = 0; i < e->n; i++)
	//	fprintf(stderr, "a[%d] = %g %g\n", i, a[i][0], a[i][1]);


	// for each of the n points, find its N nearest neighbors (in angle)
	// the zeroth point is itself
	//
	int t[e->n][e->N]; // both indices wrt the sorted order
	for (int i = 0; i < e->n; i++)
	{
		//if(i==e->i)fprintf(stderr, "i=%d\n", i);
		int n = 0;      // total number of points so far
		int p = i + 1;  // current candidate in the PLUS direction
		int m = i - 1;  // current candidate in the MINUS direction
		t[i][n++] = i;  // add i as the first point
		while (n < e->N && n < e->n)
		{
			//if(i==e->i)fprintf(stderr, "\tn=%d p=%d m=%d ",n,p,m);
			if (m < 0) m += e->n;
			if (p >= e->n) p -= e->n;
			//if (i==e->i)fprintf(stderr, "\t(p=%d m=%d) ", p, m);
			assert(m >= 0); assert(p >= 0);
			assert(m < e->n); assert(p < e->n);
			float dp = angle_dist(a[i][0], a[p][0]);
			float dm = angle_dist(a[i][0], a[m][0]);
			//if(i==e->i)fprintf(stderr, "\t{ai=%g ap=%g am=%g}",
			//		a[i][0],
			//		a[p][0],
			//		a[m][0]
			//		);
			//if(i==e->i)fprintf(stderr, "\t(dp=%g dm=%g)\n", dp, dm);
			if (dp <= dm)
				t[i][n++] = p++;
			else
				t[i][n++] = m--;
		}
		//if(i==e->i) {
		//	for (int j = 0; j < e->N; j++)
		//	{
		//		int tij = t[i][j];
		//		fprintf(stderr, "\tt[%d][%d] = %d\t"
		//				"a[%d] = %g %g\n",
		//			       	i, j, tij,
		//				tij,
		//				a[tij][0],
		//				a[tij][1]);
		//	}
		//}
		if (a[i][1]==e->i) {
			int ii = a[i][1];
			assert(ii >= 0);
			assert(ii < e->n);
			for (int j = 0; j < e->N; j++)
			if (j < 100)
			{
				int tij = t[i][j];
				int k = a[tij][1];
				e->Ni[j][0] = e->x[2*k+0];
				e->Ni[j][1] = e->x[2*k+1];
			}
		}
	}
	//for (int i = 0; i < e->n; i++)
	//for (int j = 0; j < e->N; j++)
	//{
	//	int x = t[i][j];
	//	assert(x >= 0);
	//	assert(x < e->n);
	//}

	// compute the energy of each point
	float E[e->n];
	for (int i = 0; i < e->n; i++)
	{
		int im = 0; // neighbor index of closest neighbor of i
		int iM = 0; // neighbor index of farthest neighbor of i
		for (int j = 0; j < e->N; j++)
		{
			// global indices
			int gj = a[t[i][j ]][1];
			int gm = a[t[i][im]][1];
			int gM = a[t[i][iM]][1];
			float dj = point_dist(e->c, e->x + 2*gj);
			float dm = point_dist(e->c, e->x + 2*gm);
			float dM = point_dist(e->c, e->x + 2*gM);
			if (dj < dm) im = j;
			if (dj > dM) iM = j;
			//if(a[i][1]==e->i)
			//fprintf(stderr, "j=%d dj=%g dm=%g dM=%g im=%d iM=%d\n",
			//		j, dj, dm, dM, im, iM);
		}
		float di0 = point_dist(e->c, e->x + 2*(int)a[ i      ][1]);
		float di = point_dist(e->c, e->x + 2*(int)a[t[i][0 ]][1]);
		float dm = point_dist(e->c, e->x + 2*(int)a[t[i][im]][1]);
		float dM = point_dist(e->c, e->x + 2*(int)a[t[i][iM]][1]);
		assert(di == di0);
		//if(a[i][1]==e->i)
		//fprintf(stderr, "i=%d di dm dM = %g %g %g\n", i, di, dm, dM);
		assert(di <= dM);
		assert(dm <= di);
		float x = (di - dm) / (dM - dm);
		E[i] = exp(-x*x);
		//fprintf(stderr, "E[%d] = %g\n", i, E[i]);
	}

	// fill-in the passing points in angular order
	e->m = 0;
	for (int i = 0; i < e->n; i++)
	if (E[i] > e->alpha)
	{
		int j = a[i][1]; // original index
		e->P[2*e->m + 0] = e->x[2*j+0];
		e->P[2*e->m + 1] = e->x[2*j+1];
		e->m += 1;
	}
	fprintf(stderr, "m=%d/%d\n", e->m, e->n);
}


// SECTION 5. Coordinate Conversions                                        {{{1

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

// draw a circle given center and radius
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
		f->rgb[3*idx+0] = 180;
		f->rgb[3*idx+1] = 180;
		f->rgb[3*idx+2] = 180;
	}
}

// auxiliary function for drawing a black pixel
static void plot_pixel_black(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
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

// function to draw a black segment
static void plot_segment_black(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_black, f);
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

// draw the view point
static void draw_view_center(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	float P[2];
	map_view_to_window(e, P, e->c);

	// grey circle
	float side = DISK_RADIUS;
	uint8_t green[3] = {0, 127, 0};
	splat_disk(f->rgb, f->w, f->h, P, DISK_RADIUS, green);
	//for (int j = -side-1 ; j <= side+1; j++)
	//for (int i = -side-1 ; i <= side+1; i++)
	//if (hypot(i, j) < side)
	//{
	//	int ii = P[0] + i;
	//	int jj = P[1] + j;
	//	if (insideP(f, ii, jj))
	//		for (int c = 0; c < 3; c++)
	//			f->rgb[3*(f->w*jj+ii)+c] = 127;
	//}

	// central green dot
	int ii = P[0];
	int jj = P[1];
	if (insideP(f, ii, jj))
		f->rgb[3*(f->w*jj+ii)+1]=255;
}

static void draw_x_points_in_red(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	uint8_t red[3] = {255, 0, 0};
	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		map_view_to_window(e, P, e->x + 2*i);
		splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, red);
		//plot_pixel_red(P[0], P[1], f);
		//plot_pixel_red(P[0]+1, P[1], f);
		//plot_pixel_red(P[0]-1, P[1], f);
		//plot_pixel_red(P[0], P[1]+1, f);
		//plot_pixel_red(P[0], P[1]-1, f);
	}
}

static void draw_y_points_in_gray(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	uint8_t gray[3] = {120, 120, 120};
	for (int i = 0; i < e->n; i++)
	{
		float P[2];
		map_view_to_window(e, P, e->y + 2*i);
		splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, gray);
		//plot_pixel_gray(P[0], P[1], f);
		//plot_pixel_gray(P[0]+1, P[1], f);
		//plot_pixel_gray(P[0]-1, P[1], f);
		//plot_pixel_gray(P[0], P[1]+1, f);
		//plot_pixel_gray(P[0], P[1]-1, f);
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

static void draw_c_circle_in_green(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	float P[2];
	map_view_to_window(e, P, e->c);
	plot_circle_green(f, P[0], P[1], e->r * e->scale);
}


// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state_katz(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255; // white

	draw_x_points_in_red(f);
	draw_view_center(f);
	draw_c_circle_in_green(f);

	compute_points_inversion(e);
	if (e->show_debug)
		draw_y_points_in_gray(f);

	compute_red_points_convex_hull(e);
	for (int i = 0; i < e->m - 1; i++)
	{
		float P[2], Q[2], Z[4], C[2];
		map_view_to_window(e, P, e->z + 2*i);
		map_view_to_window(e, Q, e->z + 2*i + 2);
		invert_point(Z + 0, e->c, e->r, e->z + 2*i + 0);
		invert_point(Z + 2, e->c, e->r, e->z + 2*i + 2);
		if (e->show_debug)
			plot_segment_gray(f, P[0], P[1], Q[0], Q[1]);
		map_view_to_window(e, P, Z + 0);
		map_view_to_window(e, Q, Z + 2);
		map_view_to_window(e, C, e->c);
		if (det(P, Q, C) > 0)
		{
			plot_segment_blue(f, P[0], P[1], Q[0], Q[1]);
			uint8_t blue[3] = {0, 0, 255};
			splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, blue);
		}
		else if (e->show_debug)
			plot_segment_cyan(f, P[0], P[1], Q[0], Q[1]);
	}

	// hud
	uint8_t fg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "R = %g", e->r);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, fg, NULL, 0, e->font, buf);
}

static void paint_state_biasutti(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255; // white

	// draw input data
	draw_x_points_in_red(f);
	draw_view_center(f);

	// recompute energies (always)
	biasutti(e);

	// draw visible polygon
	for (int i = 0; i < e->m; i++)
	{
		float P[2], Q[2], C[2];
		map_view_to_window(e, P, e->P + 2*i);
		map_view_to_window(e, Q, e->P + 2*((i+1)%e->m));
		map_view_to_window(e, C, e->c);
		if (det(P, Q, C) > 0)
		{
			plot_segment_blue(f, P[0], P[1], Q[0], Q[1]);
			uint8_t blue[3] = {0, 0, 255};
			splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, blue);
		}
	}

	// draw debug stuff
	if (e->show_debug)
	{
		for (int j = 0; j < e->N; j++)
		if (j < 100)
		{
			float C[2], N[2];
			map_view_to_window(e, C, e->c);
			map_view_to_window(e, N, e->Ni[j]);
			plot_segment_gray(f, C[0], C[1], N[0], N[1]);
		}
		{
			float C[2], N[2];
			map_view_to_window(e, C, e->c);
			map_view_to_window(e, N, e->x + 2*e->i);
			plot_segment_black(f, C[0], C[1], N[0], N[1]);
		}
	}

	// hud
	uint8_t fg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "N = %d\na = %g\ni = %d", e->N, e->alpha, e->i);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, fg, NULL, 0, e->font, buf);
}

static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	if (e->mode == 0)
		paint_state_katz(f);
	if (e->mode == 1)
		paint_state_biasutti(f);

	if (e->negative)
		for (int i = 0; i < 3 * f->w * f->h; i++)
			f->rgb[i] = 255 - f->rgb[i];

	f->changed = 1;
}


// SECTION 8. User-Interface Actions and Events                             {{{1

// actions to change parameters
static void shift_N(struct viewer_state *e, int d) { e->N += d; }
static void shift_i(struct viewer_state *e, int d) { e->i += d; }
static void scale_alpha(struct viewer_state *e, int d) { e->alpha += d/100.0; }

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

// action: inversion radius scale
static void change_radius(struct viewer_state *e, int x, int y, float fac)
{
	e->r *= fac;
	fprintf(stderr, "radius changed %g\n", e->r);
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
	if (k == 'd') e->show_debug = !e->show_debug;
	if (k == 'm') e->mode = !e->mode;
	if (k == 'n') e->negative = !e->negative;
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

	// N, a, i : hitboxes of font height
	// 0  1  2  3    4    5      6   7
	int Y = y / e->font->height;
	if (e->mode == 1 && k == FTR_BUTTON_DOWN && x < 10 * e->font->width)
	{
		if (Y == 0) shift_N(e, -1);
		if (Y == 1) scale_alpha(e, -1);
		if (Y == 2) shift_i(e, -1);
		f->changed = 1;
		return;
	}
	if (e->mode == 1 && k == FTR_BUTTON_UP && x < 10 * e->font->width)
	{
		if (Y == 0) shift_N(e, 1);
		if (Y == 1) scale_alpha(e, 1);
		if (Y == 2) shift_i(e, 1);
		f->changed = 1;
		return;
	}
	if (e->mode == 0 && k == FTR_BUTTON_DOWN && x < 10 * e->font->width)
	{
		if (Y == 0) change_radius(e, x, y, 1.0/RADIUS_FACTOR);
		f->changed = 1;
		return;
	}
	if (e->mode == 0 && k == FTR_BUTTON_UP && x < 10 * e->font->width)
	{
		if (Y == 0) change_radius(e, x, y, RADIUS_FACTOR);
		f->changed = 1;
		return;
	}

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

	// radius in/out (if hit), zoom in/out (if no hit)
	if (k == FTR_BUTTON_DOWN)
	{
		if (hit_point(e, x, y)<0)
			change_view_scale(e, x, y, ZOOM_FACTOR);
		else
			change_radius(e, x, y, RADIUS_FACTOR);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (hit_point(e, x, y)<0)
			change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);
		else
			change_radius(e, x, y, 1.0/RADIUS_FACTOR);
	}


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
	e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	e->x = read_ascii_floats(stdin, &e->n);
	e->y = malloc(e->n * sizeof*e->y);
	e->z = malloc(e->n * sizeof*e->y);
	e->P = malloc(e->n * sizeof*e->P);
	e->n /= 2;
	e->mode = 0; // 0=katz, 1=biasutti
	e->negative = 0;
	e->alpha = 0.95;
	e->N = 13;
	e->i = 24;
	fprintf(stderr, "read %d points from stdin\n", e->n);
	//compute_points_inversion(e);

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

int main(int c, char *v[]) { return main_pkatz(c, v); }

// vim:set foldmethod=marker:
