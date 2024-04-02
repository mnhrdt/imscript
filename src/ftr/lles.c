// lagrangian large eddy simulation (particle smoothed hydrodynamics)

#include <math.h>     // fmod, floor
#include <stdbool.h>  // bool
#include <stdio.h>    // fprintf, stdout, stderr
#include <stdlib.h>   // malloc, free, rand, RAND_MAX
#include "ftr.h"      // ftr
#include "seconds.c"  // seconds
#include "random.c"   // random_uniform, random_laplace



// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


struct les_state {
	// 1. domain data (lagrangian)
	int N;     // total number of particles
	float *P;  // particle coordinates (N*2)
	float *V;  // particle velocities  (N*2)

	// 2. domain data (eulerian)
	int w, h;     // domain dimensions
	float T0;     // temperature base (1)
	float s;      // kernel parameter
	float F;      // force multiplier
	float p;      // force exponent

	// 3. numeric parameters
	float t;      // time step

	// 4. ui
	float r;      // point radius
	struct bitmap_font font[1];
};

static void init_state(struct les_state *e, int w, int h, int n)
{
	e->w = w;
	e->h = h;
	e->N = n;
	e->P = malloc(w * n * 2 * sizeof*e->P);
	e->V = malloc(w * n * 2 * sizeof*e->V);

	e->s = 10;   // kernel parameter

	e->T0 = 0;   // base temperature
	e->t = 0.3;  // time step
	e->F = 1.0/5000; // force multiplier
	e->p = 2;

	e->r = 2.3;  // dot radius
	e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
}

static void reset_particles(struct les_state *e)
{
	// uniformly distributed particles
	for (int i = 0; i < e->N; i++) e->P[2*i + 0] = e->w * random_uniform();
	for (int i = 0; i < e->N; i++)
		e->P[2*i + 1] = e->h*1.0/3 + e->h * random_uniform()/3;

	// initial horizontal speed
	for (int i = 0; i < e->N; i++) e->V[2*i + 0] = 1;
	for (int i = 0; i < e->N; i++) e->V[2*i + 1] = 0;
}

static float fmod2(float x, float m)
{
	return x - m * floor(x / m);
}

static void move_particles(struct les_state *e)
{
	// compute force field at each particle position
	float F[e->N][2];
	for (int i = 0; i < e->N; i++)
	{
		F[i][0] = F[i][1] = 0;
		float *P = e->P + 2*i; // particle position
		float *V = e->V + 2*i; // particle velocity
		for (int j = 0; j < e->N; j++)
		if (j != i)
		{
			float *Q = e->P + 2*j;
			float D[2] = {Q[0] - P[0], Q[1] - P[1]};
			float d = fmin(1, hypot(D[0], D[1]));
			if (d < e->s)
			{
				F[i][0] -= D[0] / pow(d,e->p);
				F[i][1] -= D[1] / pow(d,e->p);
			}
		}
	}

	// apply forces to all particles
	for (int i = 0; i < e->N; i++)
	{
		float *P = e->P + 2*i; // particle position
		float *V = e->V + 2*i; // particle velocity
		V[0] += e->t * e->F * F[i][0];
		V[1] += e->t * e->F * F[i][1];
		P[0] += e->t * V[0] + e->T0*(random_normal());
		P[1] += e->t * V[1] + e->T0*(random_normal());
		P[0] = fmod2(P[0], e->w);
		P[1] = fmod2(P[1], e->h);
	}
}




static bool insideP(int w, int h, int x, int y)
{
	return  x >= 0  &&  y >= 0  &&  x < w  &&  y < h;
}

static void splat_disk(uint8_t *rgb, int w, int h, float p[2], float r,
		uint8_t color[3])
{
	for (int j = -r-1 ; j <= r+1; j++)
	for (int i = -r-1 ; i <= r+1; i++)
	if (hypot(i, j) < r)
	{
		int ii = p[0] + i;
		int jj = p[1] + j;
		float R = hypot(ii - p[0], jj - p[1]);
		if (R >= r) continue;
		if (insideP(w, h, ii, jj))
		{
			float a = pow(R/r, 2);
			for (int k = 0; k < 3; k++)
				//rgb[3*(w*jj+ii)+k] = color[k];
				//rgb[3*(w*jj+ii)+k] = a*255 + (1-a)*color[k];
				rgb[3*(w*jj+ii)+k] = a*rgb[3*(w*jj+ii)+k] + (1-a)*color[k];
		}
	}
}

// CALLBACK : idle
static void step(struct FTR *f, int x, int y, int k, int m)
{
	struct les_state *e = f->userdata;
	move_particles(e);
	//evolve_fields(e);

	// black background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 0;


	// put each particle as a white dot
	uint8_t red[3] = {255, 0, 0};
	uint8_t green[3] = {0, 255, 0};
	uint8_t white[3] = {255, 255, 255};
	uint8_t black[3] = {0, 0, 0};
	for (int i = 0; i < e->N; i++)
	{
		//int x = lrint(e->p[2*i+0]);
		//int y = lrint(e->p[2*i+1]);
		float *P = e->P + 2 * i;
		float r = e->r;
		uint8_t *color = white;
		if (i == 0) color = red;
		if (i == 1) color = green;
		splat_disk(f->rgb, f->w, f->h, P, r, color);
	}

	// hud
	uint8_t fg[3] = {0, 255, 0};
	//uint8_t bg[3] = {0, 0, 0};
	char buf[0x200];
	snprintf(buf, 0x200, "t = %g\nT0 = %g\ns = %g\n"
			"F = %g\np = %g\nr = %g\n",
			e->t, e->T0, e->s, e->F, e->p, e->r);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, fg, NULL, 0, e->font, buf);

	// invert palette
	//for (int i = 0; i < f->w * f->h * 3; i++)
	//	f->rgb[i] = 255 - f->rgb[i];

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
static void scale_pointradius(struct les_state *e, float f) { e->r *= f; }
static void scale_force(struct les_state *e, float f) { e->F *= f; }
static void scale_power(struct les_state *e, float f) { e->p *= f; }
static void shift_base_temp(struct les_state *e, float f)
{
	e->T0 += f;
	if (e->T0 < 0) e->T0 = 0;
}
static void scale_sigma(struct les_state *e, float f) { e->s *= f; }

// CALLBACK : mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct les_state *e = f->userdata;

	// t, T0, s, F, p, r,   hitboxes of font height
	// 0  1   2  3  4  5
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) scale_timestep(e, 1/1.2);
		if (Y == 1) shift_base_temp(e, -0.01);
		if (Y == 2) scale_sigma(e, 1/1.2);
		if (Y == 3) scale_force(e, 1/1.1);
		if (Y == 4) scale_power(e, 1/1.1);
		if (Y == 5) scale_pointradius(e, 1/1.1);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) scale_timestep(e, 1.2);
		if (Y == 1) shift_base_temp(e, 0.01);
		if (Y == 2) scale_sigma(e, 1.2);
		if (Y == 3) scale_force(e, 1.1);
		if (Y == 4) scale_power(e, 1.1);
		if (Y == 5) scale_pointradius(e, 1.1);
	}

	if (k == FTR_BUTTON_RIGHT) reset_particles(e);

	f->changed = 1;
}

// display another animation
int main_lles(int c, char *v[])
{
	struct les_state e[1];
	init_state(e, 800, 600, 1000);
	reset_particles(e);

	struct FTR f = ftr_new_window(e->w, e->h);
	f.userdata = e;
	ftr_set_handler(&f, "idle", step);
	ftr_set_handler(&f, "key", key);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "resize", resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

int main(int c, char *v[])
{
	return main_lles(c, v);
}
