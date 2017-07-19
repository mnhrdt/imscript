#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"


// radius of the disks that are displayed around control points
#define DISK_RADIUS 7

struct image_rgba {
	int w, h;
	uint8_t *rgba;
};

static bool insideP(struct image_rgba *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}

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
	struct image_rgba *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgba[4*idx+0] = 255;
		f->rgba[4*idx+1] = 0;
		f->rgba[4*idx+2] = 0;
		f->rgba[4*idx+3] = 255;
	}
}

// function to draw a red segment
static void plot_segment_red(struct image_rgba *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

static void add_redline(uint8_t *x, int w, int h, int from[2], int to[2])
{
	struct image_rgba f[1] = {{.w = w, .h = h, .rgba = x}};
	plot_segment_red(f, from[0], from[1], to[0], to[1]);
}

//
//// draw four red segments connecting the control points
//static void draw_four_red_segments(struct FTR *f)
//{
//	struct viewer_state *e = f->userdata;
//
//	for (int p = 0; p < 4; p++)
//	{
//		double x0[2], xf[2];
//		map_view_to_window(e, x0, e->c[p]);
//		map_view_to_window(e, xf, e->c[(p+1)%4]);
//		plot_segment_red(f, x0[0], x0[1], xf[0], xf[1]);
//	}
//}
//
//// draw disks around the control points
//static void draw_four_control_points(struct FTR *f)
//{
//	struct viewer_state *e = f->userdata;
//
//	for (int p = 0; p < 4; p++)
//	{
//		if (p == 2 && e->restrict_to_affine) continue;
//		double P[2];
//		map_view_to_window(e, P, e->c[p]);
//
//		// grey circle
//		int side = DISK_RADIUS;
//		for (int j = -side-1 ; j <= side+1; j++)
//		for (int i = -side-1 ; i <= side+1; i++)
//		if (hypot(i, j) < side)
//		{
//			int ii = P[0] + i;
//			int jj = P[1] + j;
//			if (insideP(f, ii, jj))
//				for (int c = 0; c < 3; c++)
//					f->rgb[3*(f->w*jj+ii)+c] = 127;
//		}
//
//		// central green dot
//		int ii = P[0];
//		int jj = P[1];
//		if (insideP(f, ii, jj))
//			f->rgb[3*(f->w*jj+ii)+1]=255;
//	}
//}

static void get_colour(uint8_t rgb[4], char *colour)
{
	struct caca { char *name; uint8_t r; uint8_t g; uint8_t b; } t[] = {
		{"red", 255, 0, 0},
		{"green", 0, 128, 0},
		{"gray", 127, 127, 127},
		{"blue", 0, 0, 255},
		{"white", 255, 255, 255},
		{"black", 0, 0, 0},
		{NULL, 0, 0, 0}
	};
	rgb[0] = rgb[1] = rgb[2] = 0;
	for (int i = 0; t[i].name; i++)
		if (0 == strcmp(t[i].name, colour))
		{
			fprintf(stderr, "GOT COLOR \"%s\"\n", colour);
			rgb[0] = t[i].r;
			rgb[1] = t[i].g;
			rgb[2] = t[i].b;
		}
}

static void add_dot(uint8_t *x, int w, int h, int p, int q, char *colour)
{
	struct image_rgba f[1] = {{.w = w, .h = h, .rgba = x}};
	uint8_t k[4];
	get_colour(k, colour);
	fprintf(stderr, "k = %d %d %d\n", k[0], k[1], k[2]);

	// grey circle
	int side = DISK_RADIUS;
	for (int j = -side-1 ; j <= side+1; j++)
	for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = p + i;
			int jj = q + j;
			if (insideP(f, ii, jj)) {
				for (int c = 0; c < 3; c++)
					f->rgba[4*(f->w*jj+ii)+c] = k[c];
				f->rgba[4*(f->w*jj+ii)+3] = 255;
			}
		}

	// central green dot
	int ii = p;
	int jj = q;
	if (insideP(f, ii, jj)) {
		f->rgba[4*(f->w*jj+ii)+1] = 255;
		f->rgba[4*(f->w*jj+ii)+3] = 255;
	}

}

// main function
#include "pickopt.c"
int main(int c, char *v[])
{
	bool draw_polygon = pick_option(&c, &v, "p", NULL);
	if (c < 4) {
		fprintf(stderr, "usage:\n\t%s w h o.png {xi yi colori}+\n", *v);
		//                          0 1 2 3      4  5  6...
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	char *filename_out = v[3];
	int n = (c - 4)/3;

	uint8_t *x = malloc(w * h * 4);

	for (int i = 0; i < w * h * 4; i++)
		x[i] = 50;

	if (draw_polygon) for (int i = 0; i < n; i++)
	{
		int ii = i;
		int jj = (i+1)%n;
		int from[2] = {atoi(v[4+3*ii+0]),atoi(v[4+3*ii+1])};
		int to[2]   = {atoi(v[4+3*jj+0]),atoi(v[4+3*jj+1])};
		add_redline(x, w, h, from, to);
	}
	for (int i = 0; i < n; i++)
	{
		int p = atoi(v[4+3*i+0]);
		int q = atoi(v[4+3*i+1]);
		char *colour = v[4+3*i+2];
		add_dot(x, w, h, p, q, colour);
	}

	iio_write_image_uint8_vec(filename_out, x, w, h, 4);

	return 0;
}
