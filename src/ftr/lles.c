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

	e->r = 2.3;  // dot radius
	e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
}

static void fill_in_fields(struct les_state *e)
{
	// uniformly distributed particles
	for (int i = 0; i < e->N; i++) e->P[2*i + 0] = e->w * random_uniform();
	for (int i = 0; i < e->N; i++) e->P[2*i + 1] = e->h * random_uniform();

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
	for (int i = 0; i < e->N; i++)
	{
		float *P = e->P + 2*i; // particle position
		float *V = e->V + 2*i; // particle velocity
		float F[2] = {0, 0};   // force field at this particle
		for (int j = 0; j < e->N; j++)
		if (j != i)
		{
			float *Q = e->P + 2*j;
			float D[2] = {Q[0] - P[0], Q[1] - P[1]};
			float d = hypot(D[0], D[1]);
			if (d < e->s)
			{
				F[0] -= D[0] / fmin(0.1,(d*d*d));
				F[1] -= D[1] / fmin(0.1,(d*d*d));
			}
		}
		V[0] += F[0]/5000;
		V[1] += F[1]/5000;
		P[0] += e->t * V[0] + e->T0*(random_normal());
		P[1] += e->t * V[1] + e->T0*(random_normal());
		P[0] = fmod2(P[0], e->w);
		P[1] = fmod2(P[1], e->h);
	}
}


//typedef float (*getpixel_operator)(float*,int,int,int,int);
//
//static float getpixel_abort(float *x, int w, int h, int i, int j)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h)
//		abort();
//	return x[i + j*(long)w];
//}
//
//static float getpixel_0(float *x, int w, int h, int i, int j)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h)
//		return 0;
//	return x[i + j*(long)w];
//}
//
//static float getpixel_1(float *x, int w, int h, int i, int j)
//{
//	if (i < 0) i = 0;
//	if (j < 0) j = 0;
//	if (i >= w) i = w-1;
//	if (j >= h) j = h-1;
//	return x[i+j*(long)w];
//}
//
//// like n%p, but works for all numbers
//static int gmod(int x, int m)
//{
//	int r = x % m;
//	return r < 0 ? r + m : r;
//}
//
//static int positive_reflex(int n, int p)
//{
//	int r = gmod(n, 2*p);
//	if (r == p) r -= 1;
//	if (r > p)
//		r = 2*p - r;
//	return r;
//}
//
//// extrapolate by reflection
//static float getpixel_2(float *x, int w, int h, int i, int j)
//{
//	i = positive_reflex(i, w);
//	j = positive_reflex(j, h);
//	return getpixel_abort(x, w, h, i, j);
//}
//
//static getpixel_operator getpixel = getpixel_1;
//
//static float diff_x(float *x, int w, int h, int i, int j)
//{
//	float p = getpixel(x, w, h, i+1, j);
//	float q = getpixel(x, w, h, i-1, j);
//	return (p - q)/2;
//}
//
//static float diff_y(float *x, int w, int h, int i, int j)
//{
//	float p = getpixel(x, w, h, i, j+1);
//	float q = getpixel(x, w, h, i, j-1);
//	return (p - q)/2;
//}
//
//static float laplacian(float *x, int w, int h, int i, int j)
//{
//	float a = getpixel(x, w, h, i+1, j  );
//	float b = getpixel(x, w, h, i  , j+1);
//	float c = getpixel(x, w, h, i-1, j  );
//	float d = getpixel(x, w, h, i  , j-1);
//	float p = getpixel(x, w, h, i, j);
//	return a + b + c + d - 4 * p;
//}

//static void evolve_fields(struct les_state *e)
//{
//#define FORI(n) for (int i = 0; i < n; i++)
//#define FORJ(n) for (int j = 0; j < n; j++)
//#define FORK(n) for (int j = 0; j < n; j++)
//
//	// update velocity and concentration fields according to NS equations
//	float *U = malloc(e->w * e->h * sizeof*U);
//	float *V = malloc(e->w * e->h * sizeof*V);
//	float *Q = malloc(e->w * e->h * sizeof*Q);
//	FORJ(e->h) FORI(e->w)
//	{
//		float u =   getpixel(e->u, e->w, e->h, i, j);
//		float v =   getpixel(e->v, e->w, e->h, i, j);
//		float ux =    diff_x(e->u, e->w, e->h, i, j);
//		float vx =    diff_x(e->v, e->w, e->h, i, j);
//		float Px =    diff_x(e->P, e->w, e->h, i, j);
//		float uy =    diff_y(e->u, e->w, e->h, i, j);
//		float vy =    diff_y(e->v, e->w, e->h, i, j);
//		float Py =    diff_y(e->P, e->w, e->h, i, j);
//		float uL = laplacian(e->u, e->w, e->h, i, j);
//		float vL = laplacian(e->v, e->w, e->h, i, j);
//		float q =   getpixel(e->Q, e->w, e->h, i, j); // concentration
//		float qx =    diff_x(e->Q, e->w, e->h, i, j);
//		float qy =    diff_y(e->Q, e->w, e->h, i, j);
//		float qL = laplacian(e->Q, e->w, e->h, i, j);
//		float nu = e->c1 + e->c2 * (uy + vx); // viscosity
//		U[i+j*e->w] = u + e->tau * (nu*uL - Px/e->rho - u*ux - v*uy);
//		V[i+j*e->w] = v + e->tau * (nu*vL - Py/e->rho - u*vx - v*vy);
//		Q[i+j*e->w] = q + e->tau * (e->alpha * qL - u*qx - v*qy);
//	}
//	FORI(e->w * e->h) e->u[i] = U[i];
//	FORI(e->w * e->h) e->v[i] = V[i];
//	FORI(e->w * e->h) e->Q[i] = Q[i];
//	free(Q); free(V); free(U);
//
//	// update pressure field (a few gauss-seidel iterations)
//	FORK(3)
//	FORJ(e->h) FORI(e->w)
//	{
//		float ux = diff_x(e->u, e->w, e->h, i, j);
//		float vx = diff_x(e->v, e->w, e->h, i, j);
//		float uy = diff_y(e->u, e->w, e->h, i, j);
//		float vy = diff_y(e->v, e->w, e->h, i, j);
//		float f  = - e->rho * (ux*ux + vy*vy + 2*uy*vx);
//		float p = getpixel(e->P, e->w, e->h, i, j);
//		float pL = laplacian(e->P, e->w, e->h, i, j);
//		e->P[i+j*e->w] += e->tau * (pL - f);
//	}
//	FORI(e->w * e->h)
//		e->P[i] = fmax(e->P[i], 0);
//}



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
	snprintf(buf, 0x200, "t = %g\nT0 = %g\ns = %g\nr = %g\n",
			e->t, e->T0, e->s, e->r);
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

	// t, T0, s, r  hitboxes of font height
	// 0  1   2  3
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) scale_timestep(e, 1/1.2);
		if (Y == 1) shift_base_temp(e, -0.01);
		if (Y == 2) scale_sigma(e, 1/1.2);
		if (Y == 3) scale_pointradius(e, 1/1.1);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) scale_timestep(e, 1.2);
		if (Y == 1) shift_base_temp(e, 0.01);
		if (Y == 2) scale_sigma(e, 1.2);
		if (Y == 3) scale_pointradius(e, 1.1);
	}

	f->changed = 1;
}

// display another animation
int main_lles(int c, char *v[])
{
	struct les_state e[1];
	init_state(e, 800, 600, 1000);
	fill_in_fields(e);

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
