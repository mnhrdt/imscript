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
	int w, h;  // domain dimensions
	float T0;  // temperature base (1)
	float s;   // kernel parameter
	float F;   // force multiplier
	float p;   // force exponent
	float R;   // momentum transfer coefficient

	float H;   // relative height of the input pipe
	float V0;  // initial speed at the input pipe

	int o;     // presence of obstacle
	float Oh;  // relative height of the obstacle
	float Ox;  // relative x-pos of the obstacle
	float Oy;  // relative y-pos of the obstacle


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

	e->s = 6.9;   // kernel parameter

	e->T0 = 0;   // base temperature
	e->t = 2.4;  // time step
	e->F = 0.112; // force multiplier
	e->p = 1;     // force decay
	e->R = 1.1;   // momentum drag

	e->H = 0.45;
	e->V0 = 4;

	e->o = 1;
	e->Oh = 0.07;
	e->Ox = 0.35;
	e->Oy = 0.5;

	e->r = 3.9;  // dot radius
	//e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	e->font[0] = reformat_font(*xfont_9x18B, UNPACKED);
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

// auxiliary function to compute the orientation of a triangle
// 1 : positive (anticlockwise)
// 0 : the three points are aligned
// -1 : negative (clockwise)
static int ori(float A[2], float B[2], float C[2])
{
	float U[2] = {B[0] - A[0], B[1] - A[1]};
	float V[2] = {C[0] - A[0], C[1] - A[1]};
	float d = U[0] * V[1] - U[1] * V[0];
	if (d > 0) return 1;
	if (d < 0) return -1;
	return 0;
}

// whether the segment joining two points crosses the obstacle
static bool traversed_obstacleP(struct les_state *e, float *P, float *Q)
{
	float A[2] = {e->w * e->Ox, e->h * (e->Oy - e->Oh) };
	float B[2] = {e->w * e->Ox, e->h * (e->Oy + e->Oh) };

	int ABP = ori(A, B, P);
	int ABQ = ori(A, B, Q);
	int PQA = ori(P, Q, A);
	int PQB = ori(P, Q, B);

	return ABP != ABQ && PQA != PQB;
}

// update Q as the rebound of P-->Q against the obstacle
static void reflect_against_obstacle(struct les_state *e, float *Q, float *V)
{
	// we assume that the obstacle is a vertical segment
	// thus, only the horizontal position of Q needs to be changed
	Q[0] = 2 * e->w * e->Ox - Q[0];
	V[0] = -V[0]/4;
	V[1] = V[1]/4;
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
			float *W = e->V + 2*j;
			float D[2] = {Q[0] - P[0], Q[1] - P[1]};
			float R[2] = {W[0] - V[0], W[1] - V[1]};
			float d = fmax(0.1, hypot(D[0], D[1]));
			if (d < e->s)
			{
				F[i][0] -= D[0] / pow(d,e->p);
				F[i][1] -= D[1] / pow(d,e->p);
				F[i][0] += e->R * R[0];
				F[i][1] += e->R * R[1];
			}
		}
	}

	// apply forces to all particles
	float P0[e->N][2] ; // (save previous position of each particle)
	for (int i = 0; i < e->N; i++)
	{
		float *P = e->P + 2*i; // particle position
		float *V = e->V + 2*i; // particle velocity
		P0[i][0] = P[0];
		P0[i][1] = P[1];
		V[0] += e->t * e->F * F[i][0];
		V[1] += e->t * e->F * F[i][1];
		P[0] += e->t * V[0] + e->T0*(random_normal());
		P[1] += e->t * V[1] + e->T0*(random_normal());
		// periodic boundary conditions
		// (disabled in favor of pipe re-entry)
		//P[0] = fmod2(P[0], e->w);
		//P[1] = fmod2(P[1], e->h);
	}

	// obstacle:
	// particles that tried to traverse the obstacle are reversed
	for (int i = 0; i < e->N; i++)
		if (e->o && traversed_obstacleP(e, P0[i], e->P + 2*i))
			reflect_against_obstacle(e, e->P + 2*i, e->V + 2*i);

	// pipe re-entry:
	// particles that fell outside the domain get put back in the pipe
	for (int i = 0; i < e->N; i++)
	{
		float *P = e->P + 2*i; // particle position
		float *V = e->V + 2*i; // particle velocity
		if (P[0]<0 || P[1]<0 || P[0]>=e->w-1 || P[1]>=e->h-1)
		{
			// relative coordinates of the pipe
			float A = 0.05; // relative x-position
			float H = e->H; // relative height
			P[0] = A * e->w;
			P[1] = ((random_uniform()-0.5)*H + 0.5) * e->h;
			//V[0] = 1;
			V[0] = e->V0 * (1 + 0.1*random_normal());
			V[1] = 0;
		}
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
	char buf[0x200] = {0};
	//snprintf(buf, 0x200, "t = %g\nT0 = %g\ns = %g\n"
	//		"F = %g\np = %g\nr = %g\nH = %g\n"
	//		"V0 = %g\nR = %g",
	snprintf(buf, 0x200,
			"t (timestep)       = %g\n"
			"T0 (temperature)   = %g\n"
			"s (kernel scale)   = %g\n"
			"F (force strength) = %g\n"
			"p (force decay)    = %g\n"
			"r (dot radius)     = %g\n"
			"H (pipe caliber)   = %g\n"
			"V0 (pipe speed)    = %g\n"
			"R (momentum drag)  = %g",
			e->t, e->T0, e->s, e->F, e->p, e->r, e->H, e->V0, e->R);
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

	struct les_state *e = f->userdata;
	if (k == 'o') e->o = !e->o;
}

static void scale_float(float *x, float f)  { *x *= f; }
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
	if (k == FTR_BUTTON_DOWN && x < 30 * e->font->width)
	{
		if (Y == 0) scale_float(&e->t, 1/1.2);
		if (Y == 1) shift_base_temp(e, -0.1);
		if (Y == 2) scale_float(&e->s, 1/1.2);
		if (Y == 3) scale_float(&e->F, 1/1.1);
		if (Y == 4) scale_float(&e->p, 1/1.1);
		if (Y == 5) scale_float(&e->r, 1/1.1);
		if (Y == 6) scale_float(&e->H, 1/1.1);
		if (Y == 7) scale_float(&e->V0,1/1.1);
		if (Y == 8) scale_float(&e->R, 1/1.2);
	}
	if (k == FTR_BUTTON_UP && x < 30 * e->font->width)
	{
		if (Y == 0) scale_float(&e->t, 1.2);
		if (Y == 1) shift_base_temp(e, 0.1);
		if (Y == 2) scale_float(&e->s, 1.2);
		if (Y == 3) scale_float(&e->F, 1.1);
		if (Y == 4) scale_float(&e->p, 1.1);
		if (Y == 5) scale_float(&e->r, 1.1);
		if (Y == 6) scale_float(&e->H, 1.1);
		if (Y == 7) scale_float(&e->V0,1.1);
		if (Y == 8) scale_float(&e->R, 1.2);
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
