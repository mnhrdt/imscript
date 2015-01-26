#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "ftr.h"
#include "iio.h"

#define DISK_RADIUS 7
#define ZOOM_FACTOR 1.43

struct viewer_state {
	// RGB image data
	float *img;
	int iw, ih;

	// geometry
	double c[4][2]; // control points in image coordinates
	double p[4][2]; // control points in window coordinates
	double offset[2], scale; // window viewport

	// dragging state
	bool dragging_point;
	bool dragging_ipoint;
	int dragged_point;

	// display options
	int interpolation_order; // 0=nearest, 1=linear, 2=bilinear, 3=bicubic
	bool tile_plane;
};

// reset the view
static void center_view(struct FTR *f)
{
	double margin = 33;
	struct viewer_state *e = f->userdata;

	e->c[0][0] = margin;
	e->c[0][1] = margin;
	e->c[1][0] = f->w - margin;
	e->c[1][1] = margin;
	e->c[2][0] = f->w - margin;
	e->c[2][1] = f->h - margin;
	e->c[3][0] = margin;
	e->c[3][1] = f->h - margin;

	e->p[0][0] = 0;
	e->p[0][1] = 0;
	e->p[1][0] = e->iw - 1;
	e->p[1][1] = 0;
	e->p[2][0] = e->iw - 1;
	e->p[2][1] = e->ih - 1;
	e->p[3][0] = 0;
	e->p[3][1] = e->ih - 1;

	e->dragging_point = false;
	e->dragging_ipoint = false;
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;

	e->interpolation_order = 0;
	e->tile_plane = false;
}


// wether a point is inside the window
static bool insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

#include "drawsegment.c"
#include "cmphomod.c"
#include "getpixel.c"
#include "bicubic.c"
#include "marching_interpolation.c"

static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

//// auxiliary function for drawing red segments
//static void plot_pixel_red_aa(int x, int y, float a, void *e)
//{
//	struct FTR *f = e;
//	if (insideP(f, x, y)) {
//		int idx = f->w * y + x;
//		f->rgb[3*idx+0] = f->rgb[3*idx+0]*(1-a) + a*255;
//		f->rgb[3*idx+1] = f->rgb[3*idx+1]*(1-a) + a*0;
//		f->rgb[3*idx+2] = f->rgb[3*idx+2]*(1-a) + a*0;
//	}
//}

// auxiliary function for drawing red segments
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

//// function to draw red segments
//static void plot_segment_red_aa(struct FTR *f,
//		double x0, double y0, double xf, double yf)
//{
//	traverse_segment_aa(x0, y0, xf, yf, plot_pixel_red_aa, f);
//}


// function to draw red segments
static void plot_segment_red(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

static double float_to_byte(double x)
{
	if (x < 0) return 0;
	if (x < 255) return x;
	return 255;
}

// change from view coordinates to window coordinates
static void map_view_to_window(struct FTR *f, double y[2], double x[2])
{
	struct viewer_state *e = f->userdata;
	for (int k = 0; k < 2; k++)
		y[k] = e->offset[k] + e->scale * x[k];
}

// change from window coordinates to view coordinates
static void map_window_to_view(struct FTR *f, double y[2], double x[2])
{
	struct viewer_state *e = f->userdata;
	for (int k = 0; k < 2; k++)
		y[k] = ( x[k] - e->offset[k] ) / e->scale;
}

// change from window coordinates to image coordinates
static void map_window_to_image(struct FTR *f, double y[2], double x[2])
{
	struct viewer_state *e = f->userdata;

	double H[3][3], C[4][2];
	for (int p = 0; p < 4; p++)
		map_view_to_window(f, C[p], e->c[p]);
	homography_from_4corresp(
			C[0], C[1], C[2], C[3],
			e->p[0], e->p[1], e->p[2], e->p[3],
			H);
	apply_homography(y, H, x);
}

typedef float (*getsample_operator_t)(float*,int,int,int,int,int,int);
typedef void (*interpolator_t)(float*,float*,int,int,int,float,float,
		getsample_operator_t);

static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static void bilinear_interpolation_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q, getsample_operator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	for (int l = 0; l < pd; l++) {
		float a = pix(x, w, h, pd, ip  , iq  , l);
		float b = pix(x, w, h, pd, ip+1, iq  , l);
		float c = pix(x, w, h, pd, ip  , iq+1, l);
		float d = pix(x, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

static void linear_interpolation_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q, getsample_operator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	for (int l = 0; l < pd; l++) {
		float a = pix(x, w, h, pd, ip  , iq  , l);
		float b = pix(x, w, h, pd, ip+1, iq  , l);
		float c = pix(x, w, h, pd, ip  , iq+1, l);
		float d = pix(x, w, h, pd, ip+1, iq+1, l);
		float r = marchi(a, c, b, d, p-ip, q-iq);
		result[l] = r;
	}
}

static void nearest_neighbor_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q, getsample_operator_t pix)
{
	int ip = round(p);
	int iq = round(q);
	for (int l = 0; l < pd; l++)
		result[l] = pix(x, w, h, pd, ip, iq, l);
}

static getsample_operator_t obtain_sample_operator(struct viewer_state *e)
{
	if (e->tile_plane)
		return getsample_per;
	return getsample_nan;
}

static interpolator_t obtain_interpolator(struct viewer_state *e)
{
	if (e->dragging_point || e->dragging_ipoint) return nearest_neighbor_at;
	if (e->interpolation_order == 1) return linear_interpolation_at;
	if (e->interpolation_order == 2) return bilinear_interpolation_at;
	if (e->interpolation_order == 3) return bicubic_interpolation_boundary2;
	return nearest_neighbor_at;
}


// dump the warped and the control points
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// fill-in wite background
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255;

	// draw the warped image
	if (e->img)
	{
		getsample_operator_t    pix = obtain_sample_operator(e);
		interpolator_t       interp = obtain_interpolator(e);
		double H[3][3];
		double C[4][2], P[4][2];
		for (int p = 0; p < 4; p++)
			map_view_to_window(f, C[p], e->c[p]);
		homography_from_4corresp(
				C[0], C[1], C[2], C[3],
				e->p[0], e->p[1], e->p[2], e->p[3],
				H);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		{
			double p[2] = {i, j};
			apply_homography(p, H, p);
			double w = e->iw;
			double h = e->ih;
			double ipx[2] = {
				(p[0] - 0.5) * w / (w - 1.0),
				(p[1] - 0.5) * h / (h - 1.0)
			};
			//int ip[2] = { round(ipx[0]), round(ipx[1]) };
			float v[3];
			interp(v, e->img, e->iw, e->ih, 3, ipx[0], ipx[1], pix);
			for (int l = 0; l < 3; l++)
			{
				//v[l] = getsample_per(e->img, e->iw, e->ih,
				//		3, ip[0], ip[1],l);
				if (isfinite(v[l]))
				{
					int idx = l + 3 * (f->w * j + i);
					f->rgb[idx] = float_to_byte(v[l]);
				}
			}
		}
	}

	// draw the four red segments
	for (int p = 0; p < 4; p++)
	{
		double x0[2], xf[2];
		map_view_to_window(f, x0, e->c[p]);
		map_view_to_window(f, xf, e->c[(p+1)%4]);
		plot_segment_red(f, x0[0], x0[1], xf[0], xf[1]);
	}

	// draw the four control points
	for (int p = 0; p < 4; p++)
	{
		double P[2];
		map_view_to_window(f, P, e->c[p]);

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

	f->changed = 1;
}

// action: viewport translation
static void change_view_offset(struct FTR *f, double dx, double dy)
{
	struct viewer_state *e = f->userdata;
	e->offset[0] += dx;
	e->offset[1] += dy;
}

// action: viewport zoom
static void change_view_scale(struct FTR *f, int x, int y, double fac)
{
	struct viewer_state *e = f->userdata;
	double center[2], X[2] = {x, y};
	map_window_to_view(f, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "zoom changed %g\n", e->scale);
}


// test whether (x,y) is inside one of the four control points
static int hit_point(struct FTR *f, double x, double y)
{
	struct viewer_state *e = f->userdata;
	for (int p = 0; p < 4; p++)
	{
		double P[2];
		map_view_to_window(f, P, e->c[p]);
		if (hypot(P[0] - x, P[1] - y) < DISK_RADIUS)
			return p;
	}
	return -1;
}

// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q')
	{
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'c') center_view(f);
	if (k == 'J') change_view_offset(f, 0, -1);
	if (k == 'K') change_view_offset(f, 0, 1);
	if (k == 'H') change_view_offset(f, 1, 0);
	if (k == 'L') change_view_offset(f, -1, 0);
	if (k == 'j') change_view_offset(f, 0, -10);
	if (k == 'k') change_view_offset(f, 0, 10);
	if (k == 'h') change_view_offset(f, 10, 0);
	if (k == 'l') change_view_offset(f, -10, 0);
	if (k == FTR_KEY_DOWN ) change_view_offset(f, 0, -100);
	if (k == FTR_KEY_UP   ) change_view_offset(f, 0, 100);
	if (k == FTR_KEY_RIGHT) change_view_offset(f, -100, 0);
	if (k == FTR_KEY_LEFT) change_view_offset(f, 100, 0);
	if (k == '+') change_view_scale(f, f->w/2, f->h/2, ZOOM_FACTOR);
	if (k == '-') change_view_scale(f, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);

	if (k == 'p') e->tile_plane = !e->tile_plane;
	if (isdigit(k)) e->interpolation_order = k - '0';

	paint_state(f);
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	paint_state(f);
}

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// begin dragging a control point in the WINDOW DOMAIN
	if (k == FTR_BUTTON_LEFT)
	{
		int p = hit_point(f, x, y);
		if (p >= 0)
		{
			e->dragging_point = true;
			e->dragged_point = p;
		}
	}

	// end dragging a control point in the WINDOW DOMAIN
	if (e->dragging_point && k == -FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		e->dragging_point = false;
		double X[2] = {x, y};
		map_window_to_view(f, e->c[p], X);
		paint_state(f);
	}

	// begin dragging a control point in the IMAGE DOMAIN
	if (k == FTR_BUTTON_RIGHT)
	{
		int p = hit_point(f, x, y);
		if (p >= 0)
		{
			e->dragging_ipoint = true;
			e->dragged_point = p;
		}
	}

	// end dragging a control point in the IMAGE DOMAIN
	if (e->dragging_ipoint && k == -FTR_BUTTON_RIGHT)
	{
		int p = e->dragged_point;
		e->dragging_ipoint = false;
		double P[2], Q[2] = {x, y};
		map_window_to_image(f, P, Q);
		e->p[p][0] = P[0];
		e->p[p][1] = P[1];
		map_window_to_view(f, e->c[p], Q);
		paint_state(f);
	}

	// zoom in/out
	if (k == FTR_BUTTON_DOWN) change_view_scale(f, x, y, ZOOM_FACTOR);
	if (k == FTR_BUTTON_UP) change_view_scale(f, x, y, 1.0/ZOOM_FACTOR);
}

static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// drag WINDOW DOMAIN control point (realtime feedback)
	if (e->dragging_point && m == FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		//printf("\tdragging point %d\n", p);
		double X[2] = {x, y};
		map_window_to_view(f, e->c[p], X);
		paint_state(f);
	}

	// drag IMAGE DOMAIN control point (realtime feedback)
	if (e->dragging_ipoint && m == FTR_BUTTON_RIGHT)
	{
		int p = e->dragged_point;
		double P[2], Q[2] = {x, y};
		map_window_to_image(f, P, Q);
		e->p[p][0] = P[0];
		e->p[p][1] = P[1];
		map_window_to_view(f, e->c[p], Q);
		paint_state(f);
	}
}


// main function
int main(int c, char *v[])
{
	char *filename_in = c == 2 ? v[1] : "/tmp/lenak.png";

	struct FTR f = ftr_new_window(512,512);
	struct viewer_state e[1];
	f.userdata = e;


	int pd;
	e->img = iio_read_image_float_vec(filename_in, &e->iw, &e->ih, &pd);
	assert(pd == 3);

	center_view(&f);
	paint_state(&f);

	ftr_set_handler(&f, "key", event_key);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "resize", event_resize);

	return ftr_loop_run(&f);
}
