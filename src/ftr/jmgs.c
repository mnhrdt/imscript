// Jacobi-Maupertuis geodesics simulator
//
// potentials are of the form V(r)=r^a/a

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

// force associated to the potential
// (the gradient is computed by finite differences)
static void force(float F[2], float a, float q[2])
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

// jacobi-maupertuis conformal metric associated to the potential
static float jacobi_maupertuis(float a, float E, float r)
{
	float V = potential(a, r);
	return sqrt(2*(E - V));
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
	uint8_t black[3] = {0, 0, 0};
	uint8_t pink[3] = {255, 0, 255};
	float *P = e->x, Pwin[2];
	win_from_xy(Pwin, e, P);
	splat_disk(f->rgb, f->w, f->h, Pwin, 5.3, pink);

	// force vector
	float F[2], PF[2], PFwin[2];
	force(F, e->a, P);
	PF[0] = P[0] + e->m * F[0];
	PF[1] = P[1] + e->m * F[1];
	win_from_xy(PFwin, e, PF);
	plot_segment_pink(f, Pwin[0], Pwin[1], PFwin[0], PFwin[1]);


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
			"A (bg scale)   = %g\n",
			e->a, e->E, e->m, e->bg_A
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
	// a, E, m, A   hitboxes of font height
	// 0  1  2  3
	int Y = y / e->font->height;
	int X = x / e->font->width;
	if (k == FTR_BUTTON_DOWN && x < 30 * e->font->width)
	{
		if (Y == 0) shift_float(&e->a, -0.125);
		if (Y == 1) shift_float(&e->E, -0.125);
		if (Y == 2) shift_float(&e->m, -0.125);
		if (Y == 3) shift_float(&e->bg_A, -0.125);
	}
	if (k == FTR_BUTTON_UP && x < 30 * e->font->width)
	{
		if (Y == 0) shift_float(&e->a, 0.125);
		if (Y == 1) shift_float(&e->E, 0.125);
		if (Y == 2) shift_float(&e->m, 0.125);
		if (Y == 3) shift_float(&e->bg_A, 0.125);
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
