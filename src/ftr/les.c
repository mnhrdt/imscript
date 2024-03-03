#include <math.h>     // fmod
#include <stdio.h>    // fprintf, stdout, stderr
#include <stdlib.h>   // malloc, free, rand, RAND_MAX
#include "ftr.h"      // ftr
#include "seconds.c"  // seconds


struct les_state {
	// 1. domain data
	int w, h;  // domain dimensions
	float *u;  // velocity field (w*h*2)
	float *P;  // pressure vield (w*h*1)

	// 2. numeric parameters
	float t;  // time step

	// 2. visualization data
	int n;     // number of active particles
	float *p;  // particle coordinates
};


// piecewise affine sigmoid
static float lstep(float a, float b, float t, float x)
{
	if (x < a) return 0;
	if (x > b) return t;
	return t*(x-a)/(b-a);
}


static void draw_fire(struct FTR *f, int x, int y, int k, int m)
{
	// measure time
	static int cx = 0;
	cx += 1;
	if (0 == cx % 100) {
		static double os = 0;
		double s = seconds();
		double dif = s - os;
		double fps = 100/dif;
		fprintf(stdout, "CX = %d\t%g FPS\n", cx++, fps);
		os = s;
	}

	// build palette
	static unsigned char *pal = NULL;
	if (!pal) {
		pal = malloc(3*256);
		for (int i = 0; i < 256; i++) {
			pal[3*i+0] = lstep(0,105,255,i);
			pal[3*i+1] = lstep(60,120,255,i);
			pal[3*i+2] = lstep(150,160,255,i);
		}
	}

	int num_lines_bottom = 5;
	int num_lines_hidden = 25;

	// build buffer
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
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
	int rfac = cx < 75 ? 200 : 10;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = fmod(t[p] + rfac*(rand()/(1.0+RAND_MAX)),256);
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-1; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p+2*w+1] = (1.5*t[p-3*w] + 1.7 * t[p-2*w+1]
				+ 1.5 * t[p-4*w] + 1.9 * t[p-3*w-1]
				+ 1.0 * t[p-1*w-2]
				+1.9 * t[p-4*w+1]
			) / 9.51;
	}

	// render with palette
	for (int j = 0; j < h-num_lines_hidden; j++)
	for (int i = 0; i < w; i++)
	{
		int iidx = w*(h-j-1) + i;
		int idx = (unsigned char)(lstep(105,145,255,t[iidx]));
		int pos = w*j + i;
		f->rgb[3*pos+0] = pal[3*idx+0];
		f->rgb[3*pos+1] = pal[3*idx+1];
		f->rgb[3*pos+2] = pal[3*idx+2];
	}

	f->changed = 1;
}

static void fire_resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

static void key(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

// display another animation
int main_fire(int c, char *v[])
{
	struct FTR f = ftr_new_window(800, 600);
	ftr_set_handler(&f, "idle", draw_fire);
	ftr_set_handler(&f, "key", key);
	ftr_set_handler(&f, "button", ftr_handler_toggle_idle);
	ftr_set_handler(&f, "resize", fire_resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

int main(int c, char *v[])
{
	return main_fire(c, v);
}
