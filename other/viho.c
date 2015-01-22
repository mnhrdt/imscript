#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "ftr.h"
#include "iio.h"

#define DISK_RADIUS 7

struct {
	double c[4][2];
	bool dragging;
	int dragged;
} global_state;

static bool insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

#include "drawsegment.c"

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

static void paint_state(struct FTR *f)
{
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255;
	for (int p = 0; p < 4; p++)
	{
		double x0 = global_state.c[p+0][0];
		double y0 = global_state.c[p+0][1];
		double xf = global_state.c[(p+1)%4][0];
		double yf = global_state.c[(p+1)%4][1];
		printf("segment %g %g %g %g\n", x0, y0, xf, yf);
		plot_segment_red_aa(f, x0, y0, xf, yf);
	}
	for (int p = 0; p < 4; p++)
	{
		int side = DISK_RADIUS;
		for (int j = -side-1 ; j <= side+1; j++)
		for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = global_state.c[p][0] + i;
			int jj = global_state.c[p][1] + j;
			if (insideP(f, ii, jj))
				for (int c = 0; c < 3; c++)
					f->rgb[3*(f->w*jj+ii)+c] = 0;
		}
	}
	f->changed = 1;
}

int hit_point(double x, double y)
{
	for (int p = 0; p < 4; p++)
	{
		double px = global_state.c[p][0];
		double py = global_state.c[p][1];
		if (hypot(px - x, py - y) < DISK_RADIUS)
			return p;
	}
	return -1;
}

static void print_event_button(struct FTR *f, int k, int m, int x, int y)
{
	printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
	if (k == FTR_BUTTON_LEFT)
	{
		int p = hit_point(x, y);
		if (p >= 0)
		{
			printf("\thit point %d\n", p);
			global_state.dragging = true;
			global_state.dragged = p;
		}
	}
	if (global_state.dragging && k == -FTR_BUTTON_LEFT)
	{
		int p = global_state.dragged;
		printf("\tstopped dragging point %d\n", p);
		global_state.dragging = false;
		global_state.c[p][0] = x;
		global_state.c[p][1] = y;
		paint_state(f);
	}
}

static void print_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	//printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
	if (global_state.dragging && m == FTR_BUTTON_LEFT)
	{
		int p = global_state.dragged;
		printf("\tdragging point %d\n", p);
		global_state.c[p][0] = x;
		global_state.c[p][1] = y;
		paint_state(f);
	}
}


int main(int c, char *v[])
{
	struct FTR f = ftr_new_window(512,512);

	double margin = 33;
	global_state.c[0][0] = margin;
	global_state.c[0][1] = margin;
	global_state.c[1][0] = f.w - margin;
	global_state.c[1][1] = margin;
	global_state.c[3][0] = margin;
	global_state.c[3][1] = f.h - margin;
	global_state.c[2][0] = f.w - margin;
	global_state.c[2][1] = f.h - margin;
	global_state.dragging = false;

	paint_state(&f);

	ftr_set_handler(&f, "button", print_event_button);
	ftr_set_handler(&f, "motion", print_event_motion);
	return ftr_loop_run(&f);
}
