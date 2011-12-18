#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "marching_squares.c"
#include "marching_interpolation.c"
#include "iio.h"


#ifndef xmalloc
#define xmalloc malloc
#endif//xmalloc


// draw a segment between two points
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px < qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py < qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px+j*slope), j+py, e);
		}
	}
}

// draw a segment between two points (somewhat anti-aliased)
void traverse_segment_aa(int px, int py, int qx, int qy,
		void (*f)(int,int,float,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, 1.0, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment_aa(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px < qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i <= qx-px; i++) {
				float exact = py + i*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(i+px, whole, 1-part, e);
				f(i+px, owhole, part, e);
			}
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py < qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++) {
				float exact = px + j*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(whole, j+py, 1-part, e);
				f(owhole, j+py, part, e);
			}
		}
	}
}

static void red_pixel(int x, int y, void *ii)
{
	uint8_t (**i)[3] = ii;
	i[y][x][0] = 255;
	i[y][x][1] = 0;
	i[y][x][2] = 0;
}


static int hack_width;
static int hack_height;

static void red_pixel_aa(int x, int y, float a, void *ii)
{
	if (x < 0 || y < 0 || x > hack_width || y > hack_height) return;
	uint8_t (**i)[3] = ii;
	if (i[y][x][0] != i[y][x][1]) return;
	i[y][x][0] = i[y][x][0]*(1-a) + 255*a;
	i[y][x][1] = i[y][x][0]*(1-a);
	i[y][x][2] = i[y][x][0]*(1-a);
}
static void green_pixel_aa(int x, int y, float a, void *ii)
{
	if (x < 0 || y < 0 || x > hack_width || y > hack_height) return;
	uint8_t (**i)[3] = ii;
	if (i[y][x][0] != i[y][x][1]) return;
	i[y][x][0] = i[y][x][0]*(1-a);
	i[y][x][1] = i[y][x][0]*(1-a) + 255*a;
	i[y][x][2] = i[y][x][0]*(1-a);
}

// draw contour of an image at a single threshold level
static void overlay_level_line_in_red(uint8_t (**y)[3], int n,
					float **x, int w, int h, float t)
{
	hack_width = n*w;
	hack_height = n*h;
	int ns;
	float (*s)[2][2] = marching_squares_whole_image_float(&ns,x[0],w,h,t);
	for (int i = 0; i < ns; i++) {
		traverse_segment_aa(n*s[i][0][0], n*s[i][0][1],
				n*s[i][1][0], n*s[i][1][1],
				red_pixel_aa, y);
	}
	free(s);
}

// evaluate the interpolation at a sub-pixelic point
static float marchi_eval_at(float **x, int w, int h, float p, float q)
{
	if (p < 0 || q < 0 || p+1 >= w || q+1 >= h)
		return 0;
	int ip = floor(p);
	int iq = floor(q);
	float a = x[iq][ip];
	float b = x[iq+1][ip];
	float c = x[iq][ip+1];
	float d = x[iq+1][ip+1];
	return marchi(a, b, c, d, p - ip, q - iq);
}


// utility function: alloc a 2d matrix contiguously (wxh elements of size n)
static void *matrix_build(int w, int h, size_t n)
{
	size_t p = sizeof(void *);
	char *r = xmalloc(h*p + w*h*n);
	for (int i = 0; i < h; i++)
		*(void **)(r + i*p) = r + h*p + i*w*n;
	return r;
}


// read an environment variable to toggle backround interpolation
static int toggle_interpolated_background(void)
{
	static int reat = 0;
	static int value = 0;
	if (!reat) {
		reat = 1;
		char *s = getenv("INTERPOLATE_BACKGROUND");
		if (s)
			value = atoi(s);
	}
	return value;
}

int main(int c, char *v[])
{
	// check number of arguments
	if (c < 3) {
		fprintf(stderr, "usage:\n\t%s n t1 t2 ... <image >lines\n", *v);
		return EXIT_FAILURE;
	}

	// read command line arguments
	int n = atoi(v[1]), w, h, nt = c - 2;
	float t[nt];
	for (int i = 0; i < nt; i++)
		t[i] = 0.5+atof(v[2+i]);

	// read input image
	float **x = iio_read_image_float_matrix("-", &w, &h);

	// allocate output image
	int W = n*w - n + 1;
	int H = n*h - n + 1;
	uint8_t (**y)[3] = matrix_build(W, H, sizeof**y);

	// fill background of output image
	int nh = n/2;
	float nf = n;
	for (int j = 0; j < H; j++)
		for (int i = 0; i < W; i++) {
			uint8_t g = x[(j+nh)/n][(i+nh)/n];
			if (toggle_interpolated_background() > 0)
				g = marchi_eval_at(x, w, h, i/nf, j/nf);
			for (int k = 0; k < 3; k++)
				y[j][i][k] = g;
		}

	// draw contours
	for (int i = 0; i < nt; i++)
		overlay_level_line_in_red(y, n, x, w, h, t[i]);

	// save and exit
	iio_save_image_uint8_matrix_rgb("-", y, W, H);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
