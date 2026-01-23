// lagrangian large eddy simulation (particle smoothed hydrodynamics)

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


struct kanga_state {
	int N;       // number of points
	int k;       // number of nearest neighbors
	float x[2];  // query point
	float *P;    // list of points (2*N)
	int s;       // number of iterations per redraw

	int w, h;    // window dimensions
	float a;     // accumulator (w*h)
	float xmin, xmax, ymin, ymax;

	float r;      // point radius
	struct bitmap_font font[1];
	int pause;
};

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	int p = getpid();
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_kanga_%d_%d.png", p, c);
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(n, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}

static void init_state(struct kanga_state *e, int w, int h)
{
	e->w = w;
	e->h = h;
	e->N = 500;
	e->k = 5;
	//e->a = malloc(w * n * sizeof*e->a);
	e->P = 0;//malloc(2 * e->N * sizeof*e->P);

	e->s = 1;

	e->x[0] = 2;
	e->x[1] = 3;

	e->xmin = e->ymin = -4;
	e->xmax = e->ymax = 4;

	e->r = 3.9;  // dot radius
	//e->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	e->font[0] = reformat_font(*xfont_9x18B, UNPACKED);
	e->pause = 0;
}

static int compare_distance_to_x(const void *aa, const void *bb)
{
	static float x[2];
	if (!aa) {
		const float *b = (const float *)bb;
		x[0] = b[0];
		x[1] = b[1];
		return 0;
	}
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	float A = hypot(a[0] - x[0], a[1] - x[1]);
	float B = hypot(b[0] - x[0], b[1] - x[1]);
	return ((A > B) - (A < B));
}

static void throw_particles(struct kanga_state *e)
{
	assert(e->N > 0);
	assert(e->k > 0);
	assert(e->k <= e->N);

	for (int i = 0; i < e->N; i++)
	{
		double x1 = random_uniform();
		double x2 = random_uniform();
		e->P[2*i+0] = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
		e->P[2*i+1] = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	}

	//sort_by_distance_to_x(e->P, e->N, e->x);

	compare_distance_to_x(NULL, e->x);
	qsort(e->P, e->N, 2*sizeof*e->P, compare_distance_to_x);
}

static void win_from_xy(float *ij, struct kanga_state *e, float *xy)
{
	ij[0] = e->w * (xy[0] - e->xmin) / (e->xmax - e->xmin);
	ij[1] = e->h * (xy[1] - e->ymin) / (e->ymax - e->ymin);
}

static void xy_from_win(float *xy, struct kanga_state *e, float *ij)
{
	xy[0] = e->xmin + (ij[0] / e->w) * (e->xmax - e->xmin);
	xy[1] = e->ymin + (ij[1] / e->h) * (e->ymax - e->ymin);
}

static bool insideP(int w, int h, int x, int y)
{
	return  x >= 0  &&  y >= 0  &&  x < w  &&  y < h;
}

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

// CALLBACK : idle
static void step(struct FTR *f, int x, int y, int k, int m)
{
	static int frame_counter = 0;
	static double frame_times[10] = {0};
	frame_times[frame_counter++%10] = seconds();


	struct kanga_state *e = f->userdata;
	float P[2*e->N]; e->P = P;
	throw_particles(e);

	//move_particles(e);
	//evolve_fields(e);

	// black background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 0;

	// put each particle as a white dot
	uint8_t red[3] = {255, 0, 0};
	uint8_t green[3] = {0, 255, 0};
	uint8_t dgray[3] = {100,100,100};
	uint8_t white[3] = {255, 255, 255};
	uint8_t dblue[3] = {0, 0, 227};
	uint8_t black[3] = {0, 0, 0};
	for (int i = 0; i < e->N; i++)
	{
		//int x = lrint(e->p[2*i+0]);
		//int y = lrint(e->p[2*i+1]);
		float *P = e->P + 2 * i, Pwin[2];
		uint8_t *color = dblue;
		if (i < e->k) color = white;
		win_from_xy(Pwin, e, P);
		splat_disk(f->rgb, f->w, f->h, Pwin, e->r, color);
	}
	float xwin[2];
	win_from_xy(xwin, e, e->x);
	splat_disk(f->rgb, f->w, f->h, xwin, e->r, green);

	// hud
	uint8_t fg[3] = {0, 255, 0};
	//uint8_t bg[3] = {0, 0, 0};
	char buf[0x200] = {0};
	//snprintf(buf, 0x200, "t = %g\nT0 = %g\ns = %g\n"
	//		"F = %g\np = %g\nr = %g\nH = %g\n"
	//		"V0 = %g\nR = %g",
	snprintf(buf, 0x200,
			"N (total)      = %d\n"
			"k (nearest)    = %d\n"
			"x (query)      = %5g %5g\n"
			"r (dot radius) = %g\n",
			e->N, e->k, e->x[0], e->x[1], e->r);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, fg, NULL, 0, e->font, buf);

	double fps = get_fps(frame_times, frame_counter);
	snprintf(buf, 0x200, "FPS = %g\n", fps);
	put_string_in_rgb_image(f->rgb, f->w, f->h, 0, 0+f->h-e->font->height,
			fg, NULL, 0, e->font, buf);

	// invert palette
	//for (int i = 0; i < f->w * f->h * 3; i++)
	//	f->rgb[i] = 255 - f->rgb[i];

	f->changed = 1;
}

static void action_toggle_pause(struct FTR *f)
{
	struct kanga_state *e = f->userdata;
	e->pause = !e->pause;
	ftr_set_handler(f, "idle", e->pause ? NULL : step);
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

	struct kanga_state *e = f->userdata;
	if (k == ',') action_screenshot(f);
	if (k == 'p') action_toggle_pause(f);
}

static void scale_float(float *x, float f)  { *x *= f; }
static void shift_int(int *x, int d)  { *x += d; }

// CALLBACK : mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct kanga_state *e = f->userdata;

	// left-click : move query point
	if (k == FTR_BUTTON_LEFT)
	{
		float ij[2] = {x, y};
		xy_from_win(e->x, e, ij);
	}

	// wheel : change written parameters
	// t, T0, s, F, p, r,   hitboxes of font height
	// 0  1   2  3  4  5
	int Y = y / e->font->height;
	int X = x / e->font->width;
	if (k == FTR_BUTTON_DOWN && x < 30 * e->font->width)
	{
		if (Y == 0) shift_int(&e->N, -1);
		if (Y == 1) shift_int(&e->k, -1);
		if (Y == 2 && X < 25) scale_float(e->x+0, 1/1.2);
		if (Y == 2 && X >=25) scale_float(e->x+1, 1/1.2);
		if (Y == 3) scale_float(&e->r, 1/1.1);
	}
	if (k == FTR_BUTTON_UP && x < 30 * e->font->width)
	{
		if (Y == 0) shift_int(&e->N, 1);
		if (Y == 1) shift_int(&e->k, 1);
		if (Y == 2 && X < 25) scale_float(e->x+0, 1.2);
		if (Y == 2 && X >=25) scale_float(e->x+1, 1.2);
		if (Y == 3) scale_float(&e->r, 1.1);
	}


	f->changed = 1;
}

// display another animation
#include "pickopt.c"
int main_kanga(int c, char *v[])
{
	int w = atoi(pick_option(&c, &v, "w", "1000"));
	int h = atoi(pick_option(&c, &v, "h", "1000"));

	struct kanga_state e[1];
	init_state(e, w, h);

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
	return main_kanga(c, v);
}
