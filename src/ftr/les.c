#include <math.h>     // fmod, floor
#include <stdbool.h>  // bool
#include <stdio.h>    // fprintf, stdout, stderr
#include <stdlib.h>   // malloc, free, rand, RAND_MAX
#include "ftr.h"      // ftr
#include "seconds.c"  // seconds
#include "random.c"   // random_uniform, random_laplace


#include "bilinear_interpolation.c"

// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


struct les_state {
	// 1. domain data
	int w, h;  // domain dimensions
	float *u;  // velocity field (w*h*2)
	float *P;  // pressure vield (w*h*1)
	float *T;  // temperature field (w*h*1)
	float T0;  // temperature base (1)

	// 2. numeric parameters
	float t;  // time step

	// 3. visualization data
	int n;     // number of active particles
	float *p;  // particle coordinates (n*2)

	// 4. ui
	struct bitmap_font font[1];
};

static void init_state(struct les_state *e, int w, int h, int n)
{
	e->w = w;
	e->h = h;
	e->n = n;
	e->u = malloc(w * n * 2 * sizeof*e->u);
	e->P = malloc(w * n * 1 * sizeof*e->P);
	e->T = malloc(w * n * 1 * sizeof*e->T);
	e->p = malloc(w * n * 2 * sizeof*e->p);
	e->t = 0.3;

	e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
}

static void fill_in_fields(struct les_state *e)
{
	// horizontal velocity field
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
	for (int k = 0; k < 2; k++)
		e->u[2*(j*e->w + i) + k] = 1 - k;

	// rotating field
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
	{
		float x = (i - e->w/2.0)/200;
		float y = (j - e->h/2.0)/200;
		e->u[2*(j*e->w + i) + 0] = y;
		e->u[2*(j*e->w + i) + 1] = -x;
	}

	// constant pressure
	for (int j = 0; j < e->h; j++)
	for (int i = 0; i < e->w; i++)
		e->P[j*e->w + i] = 1;

	// zero temperature
	e->T0 = 0;
	for (int i = 0; i < e->w * e->h; i++)
		e->T[i] = 0;

	// uniformly distributed particles
	for (int i = 0; i < e->n; i++) e->p[2*i + 0] = e->w * random_uniform();
	for (int i = 0; i < e->n; i++) e->p[2*i + 1] = e->h * random_uniform();
}

static float fmod2(float x, float m)
{
	return x - m * floor(x / m);
}

static void move_particles(struct les_state *e)
{
	for (int i = 0; i < e->n; i++)
	{
		float *x = e->p + 2*i + 0; // particle position
		float u[2];
		bilinear_interpolation_vec_at(u, e->u, e->w,e->h, 2, x[0],x[1]);
		for (int k = 0; k < 2; k++)
			x[k] += u[k] * e->t + e->T0*(random_laplace());
		x[0] = fmod2(x[0], e->w);
		x[1] = fmod2(x[1], e->h);
	}
}

static void evolve_fields(struct les_state *e)
{
}



static bool insideP(int w, int h, int x, int y)
{
	return  x >= 0  &&  y >= 0  &&  x < w  &&  y < h;
}

// CALLBACK : idle
static void step(struct FTR *f, int x, int y, int k, int m)
{
	struct les_state *e = f->userdata;
	move_particles(e);
	evolve_fields(e);

	// black background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 0;

	// put each particle as a red point
	for (int i = 0; i < e->n; i++)
	{
		int x = lrint(e->p[2*i+0]);
		int y = lrint(e->p[2*i+1]);
		if (i == 0)
		{
			if (insideP(f->w, f->h, x, y))
				f->rgb[3*(y*f->w+x)+1] = 255;
		} else {
			if (insideP(f->w, f->h, x, y))
			{
				f->rgb[3*(y*f->w+x)+0] = 255;
				f->rgb[3*(y*f->w+x)+1] = 255;
				f->rgb[3*(y*f->w+x)+2] = 255;
			}
		}
	}

	// hud
	uint8_t fg[3] = {0, 255, 0};
	//uint8_t bg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "s = %g\nT0 = %g\n", e->t, e->T0);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, fg, NULL, 0, e->font, buf);
}

// piecewise affine sigmoid
static float lstep(float a, float b, float t, float x)
{
	if (x < a) return 0;
	if (x > b) return t;
	return t*(x - a)/(b - a);
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

// CALLBACK : resize
static void resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

// CALLBACK : key
static void key(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

static void scale_timestep(struct les_state *e, float f) { e->t *= f; }
static void shift_base_temp(struct les_state *e, float f)
{
	e->T0 += f;
	if (e->T0 < 0) e->T0 = 0;
}

// CALLBACK : mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct les_state *e = f->userdata;

	// s, T0 hitboxes of font height
	// 0  1
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) scale_timestep(e, 1/1.3);
		if (Y == 1) shift_base_temp(e, -0.01);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) scale_timestep(e, 1.3);
		if (Y == 1) shift_base_temp(e, 0.01);
	}

	f->changed = 1;
}

// display another animation
int main_les(int c, char *v[])
{
	struct les_state e[1];
	init_state(e, 800, 600, 1000);
	fill_in_fields(e);

	struct FTR f = ftr_new_window(e->w, e->h);
	f.userdata = e;
	ftr_set_handler(&f, "idle", step);
	ftr_set_handler(&f, "key", key);
	//ftr_set_handler(&f, "button", ftr_handler_toggle_idle);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "resize", resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

int main(int c, char *v[])
{
	return main_les(c, v);
}
