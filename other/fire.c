// c99 -Ofast fire.c ftr_x11.c -lX11
#include <stdio.h>
#include <stdlib.h>
#include "ftr.h"

static float diagonal_step(float a, float b, float x)
{
	if (x < a) return 0;
	if (x > b) return 255;
	return 255*(x-a)/(b-a);
}

static void draw_fire(struct FTR *f, int x, int y, int k, int m)
{
	static int cx = 0; cx++;

	int num_lines_bottom = 5;
	int num_lines_hidden = 25;

	// build palette
	unsigned char palette[3*256];
	for (int i = 0; i < 256; i++) {
		palette[3 * i + 0] = diagonal_step(0, 105, i);
		palette[3 * i + 1] = diagonal_step(60, 120, i);
		palette[3 * i + 2] = diagonal_step(150, 160, i);
	}

	// build buffer (and restart if the size has changed)
	static float *t = NULL;
	static int w = 0, h = 0;
	if (!f || w != f->w || h != f->h + num_lines_hidden) { 
		w = f->w;
		h = f->h + num_lines_hidden;
		if (t) free(t);
		t = malloc(w * h * sizeof*f); 
		for (int i = 0; i < w*h; i++)
			t[i] = 104;
		cx = 0;
	}

	// draw random values at the bottom
	int p = 0;
	int rfac = cx < 75 ? 200 : 10; // start faster, to fill the screen
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = (unsigned char)(t[p] + rfac*(rand()/(1.0+RAND_MAX)));
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-1; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p+2*w+1] = (1.5*t[p-3*w] + 1.7 * t[p-2*w+1] 
				+ 1.5 * t[p-4*w] + 1.9 * t[p-3*w-1]
				+ 1.0 * t[p-1*w-2]
				+1.9 * t[p-4*w+1]) / 9.51;
	}

	// render with palette
	for (int j = 0; j < h-num_lines_hidden; j++)
	for (int i = 0; i < w; i++)
	{
		unsigned char idx = diagonal_step(105, 145, t[w*(h-j-1) + i]);
		f->rgb[3 * (w*j+i) + 0] = palette[3 * idx + 0];
		f->rgb[3 * (w*j+i) + 1] = palette[3 * idx + 1];
		f->rgb[3 * (w*j+i) + 2] = palette[3 * idx + 2];
	}
}

static void fire_resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

int main(void)
{
	int w = 800;
	int h = 600;
	unsigned char *x = malloc(3*w*h);
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
	ftr_set_handler(&f, "idle", draw_fire);
	ftr_set_handler(&f, "resize", fire_resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	free(x);
	return 0;
}
