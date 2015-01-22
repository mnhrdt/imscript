#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "ftr.h"
#include "iio.h"

#define DISK_RADIUS 7

struct {
	double c[4][2];
	bool dragging_point;
	int dragged_point;
	bool dragging_bg;
	int dragged_bg_handle[2];
	double offset[2];
	float *img;
	int iw, ih;
} global_state;

static void center_global_state(struct FTR *f)
{
	double margin = 33;
	global_state.c[0][0] = margin;
	global_state.c[0][1] = margin;
	global_state.c[1][0] = f->w - margin;
	global_state.c[1][1] = margin;
	global_state.c[2][0] = f->w - margin;
	global_state.c[2][1] = f->h - margin;
	global_state.c[3][0] = margin;
	global_state.c[3][1] = f->h - margin;
	global_state.dragging_point = false;
	global_state.dragging_bg = false;
	global_state.offset[0] = 0;
	global_state.offset[1] = 0;
}


static bool insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

#include "drawsegment.c"
#include "cmphomod.c"
#include "getpixel.c"
#include "bicubic.c"

//static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
//		return NAN;
//	return x[(i+j*w)*pd + l];
//}

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

static void paint_state(struct FTR *f)
{
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255;
	for (int p = 0; p < 4; p++)
	{
		double x0 = global_state.offset[0] + global_state.c[p+0][0];
		double y0 = global_state.offset[1] + global_state.c[p+0][1];
		double xf = global_state.offset[0] + global_state.c[(p+1)%4][0];
		double yf = global_state.offset[1] + global_state.c[(p+1)%4][1];
		plot_segment_red_aa(f, x0, y0, xf, yf);
	}
	for (int p = 0; p < 4; p++)
	{
		int side = DISK_RADIUS;
		for (int j = -side-1 ; j <= side+1; j++)
		for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = global_state.offset[0] + global_state.c[p][0] + i;
			int jj = global_state.offset[1] + global_state.c[p][1] + j;
			if (insideP(f, ii, jj))
				for (int c = 0; c < 3; c++)
					f->rgb[3*(f->w*jj+ii)+c] = 0;
		}
	}
	if (global_state.img)
	{
		double H[3][3];
		double from[4][2] = {
			{0,0},
			{global_state.iw-1,0},
			{global_state.iw-1,global_state.ih-1},
			{0,global_state.ih-1},
		};
		double C[4][2];
		for (int p = 0; p < 4; p++)
		for (int k = 0; k < 2; k++)
			C[p][k] = global_state.offset[k] + global_state.c[p][k];
		homography_from_4corresp(
				C[0], C[1], C[2], C[3],
				from[0], from[1], from[2], from[3],
				H);
		for (int j = 0; j < f->h; j++)
		for (int i = 0; i < f->w; i++)
		{
			double p[2] = {i, j};
			apply_homography(p, H, p);
			double w = global_state.iw;
			double h = global_state.ih;
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
			//		getsample_per
			//		);
			//
			for (int l = 0; l < 3; l++)
			{
				v[l] = getsample_nan(global_state.img,
						global_state.iw,
						global_state.ih,
						3, ip[0], ip[1],l);
				if (isfinite(v[l]))
					f->rgb[3*(f->w*j+i)+l] = v[l];//clamp(v[l],0,255);
			}
		}
	}
	for (int p = 0; p < 4; p++)
	{
		int side = DISK_RADIUS;
		for (int j = -side-1 ; j <= side+1; j++)
		for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = global_state.offset[0]+global_state.c[p][0]+i;
			int jj = global_state.offset[1]+global_state.c[p][1]+j;
			if (insideP(f, ii, jj))
				for (int c = 0; c < 3; c++)
					f->rgb[3*(f->w*jj+ii)+c] = 0;
		}
		int ii = global_state.offset[0] + global_state.c[p][0];
		int jj = global_state.offset[1] + global_state.c[p][1];
		if (insideP(f, ii, jj))
			f->rgb[3*(f->w*jj+ii)+1]=255;
	}
	f->changed = 1;
}

int hit_point(double x, double y)
{
	for (int p = 0; p < 4; p++)
	{
		double px = global_state.offset[0] + global_state.c[p][0];
		double py = global_state.offset[1] + global_state.c[p][1];
		if (hypot(px - x, py - y) < DISK_RADIUS)
			return p;
	}
	return -1;
}

static void print_event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
	if (k == 'p')
		paint_state(f);
	if (k == 'c') {
		center_global_state(f);
		paint_state(f);
	}
}

static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	paint_state(f);
}

static void print_event_button(struct FTR *f, int k, int m, int x, int y)
{
	//printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
	if (k == FTR_BUTTON_LEFT)
	{
		int p = hit_point(x, y);
		if (p >= 0)
		{
			//printf("\thit point %d\n", p);
			global_state.dragging_point = true;
			global_state.dragged_point = p;
		} else {
			//global_state.dragging_bg = true;
			global_state.dragged_bg_handle[0] = x;
			global_state.dragged_bg_handle[1] = y;
		}
	}
	if (global_state.dragging_point && k == -FTR_BUTTON_LEFT)
	{
		int p = global_state.dragged_point;
		//printf("\tstopped dragging point %d\n", p);
		global_state.dragging_point = false;
		global_state.c[p][0] = x - global_state.offset[0];
		global_state.c[p][1] = y - global_state.offset[1];
		paint_state(f);
	}
	if (global_state.dragging_bg && k == -FTR_BUTTON_LEFT)
	{
		global_state.dragging_bg = false;
		global_state.offset[0] += x - global_state.dragged_bg_handle[0];
		global_state.offset[1] += y - global_state.dragged_bg_handle[1];
		printf("offset = %g %g\n", global_state.offset[0], global_state.offset[1]);
		paint_state(f);
	}
}

static void print_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	//printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
	if (global_state.dragging_point && m == FTR_BUTTON_LEFT)
	{
		int p = global_state.dragged_point;
		//printf("\tdragging point %d\n", p);
		global_state.c[p][0] = x - global_state.offset[0];
		global_state.c[p][1] = y - global_state.offset[1];
		paint_state(f);
	}
	if (global_state.dragging_bg && m == FTR_BUTTON_LEFT)
	{
		global_state.offset[0] = x - global_state.dragged_bg_handle[0];
		global_state.offset[1] = y - global_state.dragged_bg_handle[1];
		//printf("offset = %g %g\n", global_state.offset[0], global_state.offset[1]);
		paint_state(f);
	}
}


int main(int c, char *v[])
{
	struct FTR f = ftr_new_window(512,512);

	center_global_state(&f);
	int pd;
	global_state.img = iio_read_image_float_vec("/tmp/lenak.png",
			&global_state.iw, &global_state.ih, &pd);
	assert(pd == 3);

	paint_state(&f);

	ftr_set_handler(&f, "key", print_event_key);
	ftr_set_handler(&f, "button", print_event_button);
	ftr_set_handler(&f, "motion", print_event_motion);
	ftr_set_handler(&f, "resize", event_resize);
	return ftr_loop_run(&f);
}
