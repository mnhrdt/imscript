#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "ftr.h"
#include "iio.h"

#define DISK_RADIUS 7

struct viewer_state {
	double c[4][2];
	double p[4][2];
	bool dragging_point;
	bool dragging_ipoint;
	int dragged_point;
	bool dragging_bg;
	int dragged_bg_handle[2];
	double offset[2], scale;
	float *img;
	int iw, ih;
};

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
	e->dragging_bg = false;
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;
}


static bool insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

#include "drawsegment.c"
#include "cmphomod.c"
#include "getpixel.c"
#include "bicubic.c"

static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

static void plot_pixel_red_aa(int x, int y, float a, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = f->rgb[3*idx+0]*(1-a) + a*255;
		f->rgb[3*idx+1] = f->rgb[3*idx+1]*(1-a) + a*0;
		f->rgb[3*idx+2] = f->rgb[3*idx+2]*(1-a) + a*0;
	}
}

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

static void plot_segment_red_aa(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment_aa(x0, y0, xf, yf, plot_pixel_red_aa, f);
}

static void plot_segment_red(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

static double clamp(double x, double mi, double ma)
{
	if (x < mi) return mi;
	if (x < ma) return x;
	return ma;
}

static void map_view_to_window(struct FTR *f, double y[2], double x[2])
{
	struct viewer_state *e = f->userdata;
	for (int k = 0; k < 2; k++)
		y[k] = e->offset[k] + e->scale * x[k];
}

static void map_window_to_view(struct FTR *f, double y[2], double x[2])
{
	struct viewer_state *e = f->userdata;
	for (int k = 0; k < 2; k++)
		y[k] = ( x[k] - e->offset[k] ) / e->scale;
}

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

static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	//fprintf(stderr, "pstage C0(%g %g) P0(%g %g)\n",
	//	       	e->c[0][0], e->c[0][1],
	//	       	e->p[0][0], e->p[0][1]);
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255;
	//for (int p = 0; p < 4; p++)
	//{
	//	int side = DISK_RADIUS;
	//	for (int j = -side-1 ; j <= side+1; j++)
	//	for (int i = -side-1 ; i <= side+1; i++)
	//	if (hypot(i, j) < side)
	//	{
	//		int ii = global_state.offset[0] + global_state.c[p][0] + i;
	//		int jj = global_state.offset[1] + global_state.c[p][1] + j;
	//		if (insideP(f, ii, jj))
	//			for (int c = 0; c < 3; c++)
	//				f->rgb[3*(f->w*jj+ii)+c] = 0;
	//	}
	//}
	if (e->img)
	{
		double H[3][3];
		//double from[4][2] = {
		//	{0,0},
		//	{global_state.iw-1,0},
		//	{global_state.iw-1,global_state.ih-1},
		//	{0,global_state.ih-1},
		//};
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
			int ip[2] = { round(ipx[0]), round(ipx[1]) };
			float v[3];
			//bicubic_interpolation_boundary2(v,
			//		global_state.img,
			//		global_state.iw,
			//		global_state.ih,
			//		3, ipx[0], ipx[1],
			//		getsample_0
			//		);
			
			for (int l = 0; l < 3; l++)
			{
				v[l] = getsample_nan(e->img, e->iw, e->ih,
						3, ip[0], ip[1],l);
				if (isfinite(v[l]))
					f->rgb[3*(f->w*j+i)+l] = clamp(v[l],0,255);
			}
		}
	}
	for (int p = 0; p < 4; p++)
	{
		double x0[2], xf[2];
		map_view_to_window(f, x0, e->c[p]);
		map_view_to_window(f, xf, e->c[(p+1)%4]);
		plot_segment_red_aa(f, x0[0], x0[1], xf[0], xf[1]);
	}
	for (int p = 0; p < 4; p++)
	{
		int side = DISK_RADIUS;
		double P[2];
		map_view_to_window(f, P, e->c[p]);
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
		int ii = P[0];
		int jj = P[1];
		if (insideP(f, ii, jj))
			f->rgb[3*(f->w*jj+ii)+1]=255;
	}
	f->changed = 1;
}

static void change_view_offset(struct FTR *f, double dx, double dy)
{
	struct viewer_state *e = f->userdata;
	e->offset[0] += dx;
	e->offset[1] += dy;
	paint_state(f);
}

static void change_view_scale(struct FTR *f, int x, int y, double fac)
{
	struct viewer_state *e = f->userdata;
	double center[2], X[2] = {x, y};
	map_window_to_view(f, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "\t zoom changed %g (o=%g %g)\n", e->scale, e->offset[0], e->offset[1]);
	paint_state(f);
}


int hit_point(struct FTR *f, double x, double y)
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

static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
	if (k == 'c') {
		center_view(f);
		paint_state(f);
	}

	if (k == 'j') change_view_offset(f, 0, 1);
	if (k == 'k') change_view_offset(f, 0, -1);
	if (k == 'h') change_view_offset(f, -1, 0);
	if (k == 'l') change_view_offset(f, 1, 0);
	if (k == 'J') change_view_offset(f, 0, 10);
	if (k == 'K') change_view_offset(f, 0, -10);
	if (k == 'H') change_view_offset(f, -10, 0);
	if (k == 'L') change_view_offset(f, 10, 0);
	if (k == '+') change_view_scale(f, x, y, 1.2);
	if (k == '-') change_view_scale(f, x, y, 1.0/1.2);
}

static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	paint_state(f);
}

static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;
	//printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
	if (k == FTR_BUTTON_LEFT)
	{
		int p = hit_point(f, x, y);
		if (p >= 0)
		{
			//printf("\thit point %d\n", p);
			e->dragging_point = true;
			e->dragged_point = p;
		} //else {
		//	//global_state.dragging_bg = true;
		//	e->dragged_bg_handle[0] = x;
		//	e->dragged_bg_handle[1] = y;
		//}
	}
	if (e->dragging_point && k == -FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		//printf("\tstopped dragging point %d\n", p);
		e->dragging_point = false;
		double X[2] = {x, y};
		fprintf(stderr, "theee\n");
		map_window_to_view(f, e->c[p], X);
		//e->c[p][0] = x - e->offset[0];
		//e->c[p][1] = y - e->offset[1];
		paint_state(f);
	}
	//if (e->dragging_bg && k == -FTR_BUTTON_LEFT)
	//{
	//	e->dragging_bg = false;
	//	e->offset[0] += x - e->dragged_bg_handle[0];
	//	e->offset[1] += y - e->dragged_bg_handle[1];
	//	printf("offset = %g %g\n", e->offset[0], e->offset[1]);
	//	paint_state(f);
	//}

	if (k == FTR_BUTTON_RIGHT)
	{
		int p = hit_point(f, x, y);
		if (p >= 0)
		{
			e->dragging_ipoint = true;
			e->dragged_point = p;
		}
	}
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

	if (k == FTR_BUTTON_DOWN) change_view_scale(f, x, y, 1.2);
	if (k == FTR_BUTTON_UP) change_view_scale(f, x, y, 1/1.2);
}

static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;
	//printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
	if (e->dragging_point && m == FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		//printf("\tdragging point %d\n", p);
		double X[2] = {x, y};
		map_window_to_view(f, e->c[p], X);
		paint_state(f);
	}
	//if (e->dragging_bg && m == FTR_BUTTON_LEFT)
	//{
	//	e->offset[0] = x - e->dragged_bg_handle[0];
	//	e->offset[1] = y - e->dragged_bg_handle[1];
	//	//printf("offset = %g %g\n", e->offset[0], e->offset[1]);
	//	paint_state(f);
	//}
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
