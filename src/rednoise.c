#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

int width, height;

static void red_pixel(int i, int j, void *data)
{
	float *x = data;
	float *p = x + 3*(width*j + i);
	p[0] = 255;
	p[1] = 0;
	p[2] = 0;
}

static void red_pixel_aa(int i, int j, float a, void *data)
{
	float *x = data;
	float *p = x + 3*(width*j + i);
	p[0] = p[0]*(1-a)+255*a;
	p[1] = p[1]*(1-a);
	p[2] = p[2]*(1-a);
}

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
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
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

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

int main()
{
	int w, h;
	float *x = iio_read_image_float_rgb("-", &w, &h);

	width = w;
	height = h;
	for (int i = 0; i < 100; i++)
	{
		int p[4] = {randombounds(1, w-1), randombounds(1, h-1),
			randombounds(1, w-1), randombounds(1, h-1)};
		traverse_segment(p[0], p[1], p[2], p[3], red_pixel, x);
	}

	iio_save_image_float_vec("-", x, w, h, 3);

	return 0;
}
