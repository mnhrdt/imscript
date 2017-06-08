// c99 -Ofast foc.c ftr_x11.c -lX11
#include <stdio.h>
#include <stdlib.h>
#include "ftr.h"

static void draw_fire(struct FTR *f, int x, int y, int k, int m)
{
	int num_lines_bottom = 3;

	// build palette
	unsigned char palette[3*256];
	for (int i = 0; i < 256; i++) {
		palette[3 * i + 0] = 4 * i;
		palette[3 * i + 1] = (255-2 * i)/3;
		palette[3 * i + 2] = 3 * i;
	}

	// build buffer (if resize, restart buffer)
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
	if (!f || w != f->w || h != f->h) { 
		w = f->w;
		h = f->h;
		if (t) free(t);
		t = malloc(w * h * sizeof*t); 
		for (int i = 0; i < w*h; i++)
			t[i] = (unsigned char)rand();
	}

	// draw random values at the bottom
	int p = 0;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = (unsigned char)(t[p] + 15*(rand()/(1.0+RAND_MAX)));
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-1; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p] = (t[p-w] + 2 * t[p-2*w-1] + 2 * t[p-2*w] + 2 * t[p-2*w+1])
			* 9 * 4 / 256;
	}

	// render with palette
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = (unsigned char)(t[w*(h-j-1) + i]);
		f->rgb[3 * (w*j+i) + 0] = palette[3 * idx + 0];
		f->rgb[3 * (w*j+i) + 1] = palette[3 * idx + 1];
		f->rgb[3 * (w*j+i) + 2] = palette[3 * idx + 2];
	}
}

static void print_resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

int main(void)
{
	struct FTR f = ftr_new_window(800, 600);
	ftr_set_handler(&f, "idle", draw_fire);
	ftr_set_handler(&f, "resize", print_resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}
