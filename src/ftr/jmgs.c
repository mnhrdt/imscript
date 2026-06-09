// Jacobi-Maupertuis geodesics simulator
//
// potentials are of the form V(r)=(r^a-1)/a
// NOTE: this means that for the newtonian case a=-1, the critical energy
// corresponding to parabolic orbits is E=1, not E=0
//
// TODO:
// - draw tissot's indicators of the metric
//   Required state variables: tissot_scale, tissot_n
// - implement solver for geodesic equations (geometric leapfrog?)
// - add global toggle for mechanics/geometry modes
// DONE
// - implement symplectic euler method for trajectories q''=-grad(V)(q)
//   Required state variables: nsteps, hstep, v0_angle
// - implement symplectic leapfrog (same variables)

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


struct jmg_state {
	// physical parameters
	float a;   // potential exponent (-1=newton, 2=hooke, 0=logarithmic)
	float E;   // total energy
	float m;   // "mass" (amounts to a global rescaling)

	// domain parameters
	int w, h;  // window dimensions
	float xmin, xmax, ymin, ymax;

	// visualisation status
	float x[2];     // point of interest
	float bg_mode;  // 0=potential, 1=metric
	float bg_A;     // color scale
	int tissot_n;
	float tissot_scale;

	// trajectories and geodesics
	float j0, v0;   // initial angle and speed (cauchy data)
	int solver;     // 0,1=euler direct,symplectic 2=vel.verlet 3=yosida
	float tstep;    // timestep
	int N;          // number of steps

	// gui
	struct bitmap_font font[1];
};

// homogeneous power-law potentials of the form (r^a-1)/a
// for a=0 it is log(r)
static float potential(float a, float r)
{
	if (a == 0)
		return log(r);
	else
		return (pow(r,a) - 1)/a;
}



#define FORK(_n) for(int k=0;k<_n;k++)

// corresponding force F = -grad(V) = -r^(a-2)*q (for any a)
static void force(float F[2], float a, float q[2])
{
	float r = hypot(q[0], q[1]);
	F[0] = - pow(r, a-2) * q[0];
	F[1] = - pow(r, a-2) * q[1];
}


// force associated to the potential
// (the gradient is computed by finite differences)
static void force_finitediff(float F[2], float a, float q[2])
{
	float e = 0.00001;
	float r00 = hypot(q[0], q[1]);
	float r10 = hypot(q[0]+e, q[1]);
	float r01 = hypot(q[0], q[1]+e);
	float V00 = potential(a, r00);
	float V10 = potential(a, r10);
	float V01 = potential(a, r01);
	F[0] = -(V10 - V00)/e;
	F[1] = -(V01 - V00)/e;
}


typedef void (*vector_field)(float *, float *);

static void powerlaw_force(float *F, float *p)
{
	static float a = NAN;
	if (!F) a = *p;
	else force(F, a, p);
}


// jacobi-maupertuis conformal metric associated to the potential
// NOTE: hill region is defined by jacobi_maupertuis being NaN
static float jacobi_maupertuis(float a, float E, float r)
{
	float V = potential(a, r);
	return sqrt(2*(E - V));
}

// cutoff radius (for an embedding as a surface of revolution
static float cutoff_radius(float a, float E)
{
	if (a < -4) return 0;
	if (a == -4 && E >= 0.25) return INFINITY;
	if (a == -4 && E < 0.25) return 0;
	if (a > -4 && a < 0 && 1+a*E > 0) return pow(4*(1+a*E)/(a+4),1/a);
	if (a > -4 && a < 0 && 1+a*E <= 0) return INFINITY;
	if (a == 0) return exp(E-0.25);
	if (a > 0) return pow(4*(a*E+1)/(a+4),1/a);
	return NAN;
}



typedef void (*ode_solver)(float*,vector_field,float[2],float[2],int,float);


static void euler_direct(float *o,         // output points
		vector_field F,            // force field
		float q0[2], float v0[2],  // input data
		int N, float h)            // parameters
{
	float q[2] = { q0[0], q0[1] };
	float v[2] = { v0[0], v0[1] };
	float a[2];
	for (int i = 0; i < N; i++)
	{       // drift then kick
		F(a, q);
		FORK(2) q[k] += h * v[k];
		FORK(2) v[k] += h * a[k];
		FORK(2) o[2*i+k] = q[k];
	}
}

static void euler_symplectic(float *o,     // output points
		vector_field F,            // force field
		float q0[2], float v0[2],  // input data
		int N, float h)            // parameters
{
	float q[2] = { q0[0], q0[1] };
	float v[2] = { v0[0], v0[1] };
	float a[2];
	for (int i = 0; i < N; i++)
	{       // kick then drift
		F(a, q);
		FORK(2) v[k] += h * a[k];
		FORK(2) q[k] += h * v[k];
		FORK(2) o[2*i+k] = q[k];
	}
}


// updates q and v in-place according to velocity Verlet half-step rule
// useful for Verlet leapfrog and for Yosida solvers
static void one_verlet_step(float q[2], float v[2], vector_field F, float h)
{
	// half-kick, drift, half-kick
	float a[2];
	F(a, q);
	FORK(2) v[k] += 0.5 * h * a[k];
	FORK(2) q[k] += h * v[k];
	F(a, q);
	FORK(2) v[k] += 0.5 * h * a[k];
}


// https://en.wikipedia.org/wiki/Leapfrog_integration
static void verlet_leapfrog(float *o,     // output points
		vector_field F,            // force field
		float q0[2], float v0[2],  // input data
		int N, float h)            // parameters
{
	float q[2] = { q0[0], q0[1] };
	float v[2] = { v0[0], v0[1] };
	for (int i = 0; i < N; i++)
	{
		one_verlet_step(q, v, F, h);
		FORK(2) o[2*i+k] = q[k];
	}
}

// https://en.wikipedia.org/wiki/Leapfrog_integration#Yoshida_algorithms
static void yosida_4th_order (float *o,    // output points
		vector_field F,            // force field
		float q0[2], float v0[2],  // input data
		int N, float h)            // parameters
{
	// the famous Yosida constants
	float w0 = cbrt(2)/(cbrt(2)-2); // -1.70...
	float w1 = 1/(2-cbrt(2));       //  1.35...

	float q[2] = { q0[0], q0[1] };
	float v[2] = { v0[0], v0[1] };
	for (int i = 0; i < N; i++)
	{
		one_verlet_step(q, v, F, w1*h);
		one_verlet_step(q, v, F, w0*h);
		one_verlet_step(q, v, F, w1*h);
		FORK(2) o[2*i+k] = q[k];
	}
}


static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	int p = getpid();
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_jmgs_%d_%d.png", p, c);
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(n, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}

static void init_state(struct jmg_state *e, int w, int h)
{
	e->a = -1;
	e->E = 0;
	e->m = 1;

	e->w = w;
	e->h = h;
	e->xmin = e->ymin = -2;
	e->xmax = e->ymax = 2;

	e->bg_mode = 1;
	e->bg_A = 1;
	e->x[0] = 1;
	e->x[1] = 0;

	e->j0 = 20;
	e->v0 = 0.5;
	e->solver = 0;
	e->N = 100;
	e->tstep = 0.015625;

	e->tissot_n = 7;
	e->tissot_scale = 1;

	//e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	e->font[0] = reformat_font(*xfont_9x18B, UNPACKED);
}

static void win_from_xy(float *ij, struct jmg_state *e, float *xy)
{
	ij[0] = e->w * (xy[0] - e->xmin) / (e->xmax - e->xmin);
	ij[1] = e->h * (-xy[1] - e->ymin) / (e->ymax - e->ymin);
}

static void xy_from_win(float *xy, struct jmg_state *e, float *ij)
{
	xy[0] = e->xmin + (ij[0] / e->w) * (e->xmax - e->xmin);
	xy[1] = -(e->ymin + (ij[1] / e->h) * (e->ymax - e->ymin));
}

static bool insideP(int w, int h, int x, int y)
{
	return  x >= 0  &&  y >= 0  &&  x < w  &&  y < h;
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

// traverse a circle given center and radius
static void traverse_circle(int cx, int cy, int r,
		void (*f)(int,int,void*), void *e)
{
	int h = r / sqrt(2);
	for (int i = -h; i <= h; i++)
	{
		int s = sqrt(r*r - i*i);
		f(cx + i, cy + s, e); // upper quadrant
		f(cx + i, cy - s, e); // lower quadrant
		f(cx + s, cy + i, e); // right quadrant
		f(cx - s, cy + i, e); // left quadrant
	}
}

// auxiliary function for drawing a black pixel
static void plot_pixel_black(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f->w, f->h, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}
static void plot_pixel_pink(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f->w, f->h, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 255;
	}
}
static void plot_pixel_cyan(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f->w, f->h, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 255;
		f->rgb[3*idx+2] = 255;
	}
}



// function to draw a black segment
static void plot_segment_black(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_black, f);
}
static void plot_segment_pink(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_pink, f);
}

static void plot_circle_pink(struct FTR *f, float x, float y, float r)
{
	traverse_circle(x, y, r, plot_pixel_pink, f);
}
static void plot_circle_cyan(struct FTR *f, float x, float y, float r)
{
	traverse_circle(x, y, r, plot_pixel_cyan, f);
}

// function to draw a colored blob/disk
static void splat_disk(uint8_t *rgb, int w, int h, float p[2], float r,
		uint8_t color[3])
{
	if (r <= 1)
	{
		int ii = p[0];
		int jj = p[1];
		for (int k = 0; k < 3; k++)
			rgb[3*(w*jj+ii)+k] = color[k];
		return;
	}
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

static double get_fps(double t[10], int c)
{
	int prev = c%10 - 2;
	int curr = c%10 - 1;
	if (prev < 0) prev += 10;
	if (curr < 0) curr += 10;
	return 1.0/(t[curr] - t[prev]);
}

static void palette_jm(uint8_t *rgb, float A, float S)
{
	uint8_t black[3] = {0, 0, 0};
	uint8_t beige[3] = {0xa6, 0x7b, 0x5b};
	uint8_t white[3] = {255, 255, 255};
	if (isfinite(S))
	{
		assert(S >= 0);
		float r = S<A ? 2-S/A : A/S;
		if (r > 1) // interpolate from black to beige
			for (int k = 0; k < 3; k++)
				rgb[k] = (r-1)*white[k] + (2-r)*beige[k];
		else // interpolate from beige to white
			for (int k = 0; k < 3; k++)
				rgb[k] = r*beige[k] + (1-r)*black[k];
	} else
	{
		rgb[0] = 150;
		rgb[1] = 200;
		rgb[2] = 150;
	}
}

// CALLBACK : expose
static void event_expose(struct FTR *f, int ev_b, int ev_m, int ev_x, int ev_y)
{
	struct jmg_state *e = f->userdata;

	// gray background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 0;//127;

	// metric field
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		float p[2], wp[2] = {i, j};
		xy_from_win(p, e, wp);
		float S = jacobi_maupertuis(e->a, e->E, hypot(p[0], p[1]));
		uint8_t *rgb = f->rgb + 3*(j * f->w + i);
		palette_jm(rgb, e->bg_A, S);
	}

	// tissots
	for (int j = 0; j < e->tissot_n; j++)
	for (int i = 0; i < e->tissot_n; i++)
	{
		float pw[2] = {i*e->w*1.0/e->tissot_n, j*e->h*1.0/e->tissot_n};
		float p[2];
		xy_from_win(p, e, pw);
		float S = jacobi_maupertuis(e->a, e->E, hypot(p[0], p[1]));
		float R = e->tissot_scale/S;
		if (!isfinite(S)) continue;
		if (!isfinite(R)) continue;
		plot_circle_pink(f, pw[0], pw[1], R);
		//fprintf(stderr, "circle tissot %g %g R=%g\n", pw[0],pw[1],R);
	}

	// axes
	float ax[3][2] = {{0,0},{1,0},{0,1}};
	float wax[3][2];
	for (int i = 0; i < 3; i++)
		win_from_xy(wax[i], e, ax[i]);
	plot_segment_black(f, wax[0][0], wax[0][1], wax[1][0], wax[1][1]);
	plot_segment_black(f, wax[0][0], wax[0][1], wax[2][0], wax[2][1]);


	// point of interest
	uint8_t red[3] = {255, 0, 0};
	uint8_t green[3] = {0, 255, 0};
	uint8_t dgreen[3] = {0, 155, 0};
	uint8_t dgray[3] = {100,100,100};
	uint8_t white[3] = {255, 255, 255};
	uint8_t dblue[3] = {0, 0, 227};
	uint8_t cyan[3] = {0, 255, 255};
	uint8_t black[3] = {0, 0, 0};
	uint8_t pink[3] = {255, 0, 255};
	float *P = e->x, Pwin[2];
	win_from_xy(Pwin, e, P);
	splat_disk(f->rgb, f->w, f->h, Pwin, 5.3, pink);

	//// force vector
	//float F[2], PF[2], PFwin[2];
	//force(F, e->a, P);
	//PF[0] = P[0] + e->m * F[0];
	//PF[1] = P[1] + e->m * F[1];
	//win_from_xy(PFwin, e, PF);
	//plot_segment_pink(f, Pwin[0], Pwin[1], PFwin[0], PFwin[1]);

	// admissibility region
	float rtop = cutoff_radius(e->a, e->E);
	if (isfinite(rtop)) {
		float p[2][2] = { {0,0}, {rtop, 0} }, P[2][2];
		win_from_xy(P[0], e, p[0]);
		win_from_xy(P[1], e, p[1]);
		float R = hypot(P[1][0]-P[0][0], P[1][1]-P[0][1]);
		plot_circle_cyan(f, P[0][0], P[0][1], R);
	}

	// velocity vector
	if (true) //  if "metric" mode
	{
		float r = hypot(e->x[0], e->x[1]);
		e->v0 = sqrt(2*(e->E - potential(e->a, r))/e->m);
	}
	float θ = e->j0 * M_PI / 180;
	float V[2] = { e->v0 * cos(θ), e->v0 * sin(θ) }, PV[2], PVw[2];
	PV[0] = P[0] + V[0];
	PV[1] = P[1] + V[1];
	win_from_xy(PVw, e, PV);
	if (isfinite(e->v0))
		plot_segment_black(f, Pwin[0], Pwin[1], PVw[0], PVw[1]);
	else
		put_string_in_rgb_image(f->rgb, f->w, f->h,
			f->w/2-50, f->h/4,
			red, black, 0, e->font, " IMPOSSIBLE V0 ");

	// solver
	ode_solver solvers[4] = {
		euler_direct,
		euler_symplectic,
		verlet_leapfrog,
		yosida_4th_order
	};
	float pp[2*e->N];
	powerlaw_force(NULL, &e->a);
	if (isfinite(e->v0))
	{
		solvers[e->solver](pp, powerlaw_force, e->x, V, e->N, e->tstep);
		for (int i = 0; i < e->N; i++)
		{
			float ow[2];
			win_from_xy(ow, e, pp + 2*i);
			splat_disk(f->rgb, f->w, f->h, ow, 3.3, red);
		}

	}


	// hud
	uint8_t *hud_fg = dgreen;
	uint8_t *hud_bg = black;
	//uint8_t bg[3] = {0, 0, 0};
	char buf[0x200] = {0};
	//snprintf(buf, 0x200, "t = %g\nT0 = %g\ns = %g\n"
	//		"F = %g\np = %g\nr = %g\nH = %g\n"
	//		"V0 = %g\nR = %g",
	snprintf(buf, 0x200,
			"a (potential)  = %g\n"
			"E (energy)     = %g\n"
			"m (mass)       = %g\n"
			"A (bg scale)   = %g\n"
			"j0 (angle)     = %g\n"
			"v0 (speed)     = %g\n"
			"solver         = %d\n"
			"N (nsteps)     = %d\n"
			"h (timestep)   = %g\n",
			e->a, e->E, e->m, e->bg_A, e->j0, e->v0,
			e->solver, e->N, e->tstep
		);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, hud_fg, hud_bg, 0, e->font, buf);

	// lower hud with point data
	float r = hypot(e->x[0], e->x[1]);
	snprintf(buf, 0x200,
			"p = %g %g (%g)\n"
			"V = %g\n"
			"S = %g\n",
			e->x[0], e->x[1], r,
			potential(e->a, r),
			jacobi_maupertuis(e->a, e->E, r)
		);
	put_string_in_rgb_image(f->rgb, f->w, f->h, 0, 0+f->h-3*e->font->height,
			pink, hud_bg, 0, e->font, buf);



	//double fps = get_fps(frame_times, frame_counter);
	//snprintf(buf, 0x200, "FPS = %g\n", fps);
	//put_string_in_rgb_image(f->rgb, f->w, f->h, 0, 0+f->h-e->font->height,
	//		fg, NULL, 0, e->font, buf);

	// invert palette
	//for (int i = 0; i < f->w * f->h * 3; i++)
	//	f->rgb[i] = 255 - f->rgb[i];

	f->changed = 1;
}



// CALLBACK : resize
static void event_resize(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "resize %d %d\n", x, y);
}

// CALLBACK : key
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);

	struct jmg_state *e = f->userdata;
	if (k == ',') action_screenshot(f);
	////if (k == 'p') action_toggle_pause(f);
}

static void scale_float(float *x, float f)  { *x *= f; }
static void shift_float(float *x, float f)  { *x += f; }
static void shift_int(int *x, int d)  { *x += d; }
static void shift_angle(float *x, float d) { *x = remainder(*x + d, 360); }
static void cycle_int(int *x, int d, int m)
{
	*x += d;
	if (*x >= m) *x -= m;
	if (*x < 0) *x += m;
}
static void scale_int(int *x, float f)  { *x *= f; }

// CALLBACK : mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct jmg_state *e = f->userdata;

	// left-click : move query point
	if (k == FTR_BUTTON_LEFT)
	{
		float ij[2] = {x, y};
		xy_from_win(e->x, e, ij);
	}

	// wheel : change written parameters
	// a, E, m, A, j0, v0, solver, N, h, Tn, Ts
	// 0  1  2  3  4   5   6       7  8  9   10
	// (hitboxes of font height)
	int Y = y / e->font->height;
	int X = x / e->font->width;
	if (k == FTR_BUTTON_DOWN && x < 30 * e->font->width)
	{
		if (Y == 0) shift_float(&e->a, -0.125);
		if (Y == 1) shift_float(&e->E, -0.125);
		if (Y == 2) shift_float(&e->m, -0.125);
		if (Y == 3) shift_float(&e->bg_A, -0.125);
		if (Y == 4) shift_angle(&e->j0, -10);
		if (Y == 5) scale_float(&e->v0, 1/pow(2,0.25));
		if (Y == 6) cycle_int(&e->solver, -1, 4);
		if (Y == 7) scale_int(&e->N, 1.0/1.1);
		if (Y == 8) scale_float(&e->tstep, 1/pow(2,0.25));
		if (Y == 9) shift_int(&e->tissot_n, -1);
		if (Y ==10) scale_float(&e->tissot_scale, 1/1.1);
	}
	if (k == FTR_BUTTON_UP && x < 30 * e->font->width)
	{
		if (Y == 0) shift_float(&e->a, 0.125);
		if (Y == 1) shift_float(&e->E, 0.125);
		if (Y == 2) shift_float(&e->m, 0.125);
		if (Y == 3) shift_float(&e->bg_A, 0.125);
		if (Y == 4) shift_angle(&e->j0, 10);
		if (Y == 5) scale_float(&e->v0, pow(2,0.25));
		if (Y == 6) cycle_int(&e->solver, 1, 4);
		if (Y == 7) scale_int(&e->N, 1.1);
		if (Y == 8) scale_float(&e->tstep, pow(2,0.25));
		if (Y == 9) shift_int(&e->tissot_n, 1);
		if (Y ==10) scale_float(&e->tissot_scale, 1.1);
	}

	f->changed = 1;
}

#include "pickopt.c"
int main_jmgs(int c, char *v[])
{
	int w = atoi(pick_option(&c, &v, "w", "800"));
	int h = atoi(pick_option(&c, &v, "h", "800"));

	struct jmg_state e[1];
	init_state(e, w, h);

	struct FTR f = ftr_new_window(e->w, e->h);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "key", event_key);
	ftr_set_handler(&f, "button", event_button);
	//ftr_set_handler(&f, "resize", resize);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

int main(int c, char *v[])
{
	return main_jmgs(c, v);
}
