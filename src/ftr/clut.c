// CPU-based Glut-inspider rendered
//
// - three ingredients: frustum matrix, list of triangles, list of lights
// - no global state
// - goal: display simple polygons and surfaces interactively


#include <math.h>     // exp, pow
#include <limits.h>   // INT_MIN
#include <stdbool.h>  // bool
#include <stdint.h>   // uint8_t
#include <stdio.h>    // fprintf, stdout, stderr
#include <stdlib.h>   // malloc,free
#include <unistd.h>   // getpid
#include "ftr.h"      // ftr


// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"

struct clut_state {
	float M[4][4];  // view matrix
	float (*v)[3];  // vertex coordinates
	int nv;         // number of vertices
	int (*t)[3];    // triangle indices
	int nt;         // number of triangles

	float light[4]; // light position (projective)

	float w;        // viewport width in pixels
	float h;        // viewport height in pixels

	// controls
	float O[3];     // global offset (used to compute the view matrix)
	float a,b,c;    // euler angles (used to compute the view matrix)

	// gui
	struct bitmap_font font[1];
	int hud;
};

static void clut_emtpy(struct clut_state *e)
{
	// no mesh
	e->v = NULL;
	e->t = NULL;
	e->nv = e->nt = 0;

	// starting position = identity matrix
	e->M[0][0]=1; e->M[0][1]=0; e->M[0][2]=0; e->M[0][3]=0;
	e->M[1][0]=0; e->M[1][1]=1; e->M[1][2]=0; e->M[1][3]=0;
	e->M[2][0]=0; e->M[2][1]=0; e->M[2][2]=1; e->M[2][3]=0;
	e->M[3][0]=0; e->M[3][1]=0; e->M[3][2]=0; e->M[3][3]=1;

	// light at zenith
	e->light[0] = 0; e->light[1] = 0; e->light[2] = 0; e->light[3] = 1;

	// controls (offset + euler angles)
	e->O[0] = 0;
	e->O[1] = 0;
	e->O[2] = 0;
	e->a = 10;
	e->b = 20;
	e->c = 30;

	e->font[0] = reformat_font(*xfont_9x18B, UNPACKED);
	e->hud = 1;
}

static void clut_free(struct clut_state *e)
{
	free(e->t);
	free(e->v);
}

//static const float simplex_v[4][3] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };
//static const int simplex_T[4][3]   = { {1,2,3}, {0,1,2}, {0,2,3}, {0,3,1} };

static void clut_fill_simplex(struct clut_state *e)
{
	float V[4][3] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };
	int T[4][3]   = { {1,2,3}, {0,1,2}, {0,2,3}, {0,3,1} };
	e->nv = 4;
	e->v = malloc(e->nv*3*sizeof(float));
	e->nt = 4;
	e->t = malloc(e->nt*3*sizeof(int));
	for (int i = 0; i < e->nv; i++)
	for (int k = 0; k < 3; k++)
		e->v[i][k] = V[i][k];
	for (int i = 0; i < e->nt; i++)
	for (int k = 0; k < 3; k++)
		e->t[i][k] = T[i][k];

}

static void sqmp3(float Z[3][3], float X[3][3], float Y[3][3])
{
	float T[3][3];
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		T[i][j] = 0;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	for (int k = 0; k < 3; k++)
		T[i][j] += X[i][k] * Y[k][j];
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		Z[i][j] = T[i][j];
}

static void fill_matrix_from_controls(struct clut_state *e)
{
	float α = e->a * M_PI / 180;
	float β = e->b * M_PI / 180;
	float γ = e->c * M_PI / 180;
	float Ma[3][3] = {
		{cos(α), -sin(α), 0},
		{sin(α), cos(α), 0},
		{0, 0, 1} };
	float Mb[3][3] = {
		{cos(β), 0, sin(β)},
		{0, 1, 0},
		{-sin(β), 0, cos(β)} };
	float Mc[3][3] = {
		{1, 0, 0},
		{0, cos(γ), -sin(γ)},
		{0, sin(γ), cos(γ)} };
	float R[3][3];
	sqmp3(R, Ma, Mb);
	sqmp3(R, R, Mc);

	// TODO: build the matrix here from the viewpoint
	float (*M)[4] = e->M;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		M[i][j] = R[i][j];
	for (int i = 0; i < 3; i++)
		M[3][i] = M[i][3] = 0;
	M[3][3] = 1;

	//fprintf(stderr, "abc = %g %g %g\n", e->a, e->b, e->c);
	//fprintf(stderr, "αβγ = %g %g %g\n", α, β, γ);
	//fprintf(stderr, "M =\n");
	//for (int i = 0; i < 4; i++)
	//for (int j = 0; j < 4; j++)
	//	fprintf(stderr, "%g%c", M[i][j], j==3?'\n':' ');
}

static bool project(float *ij, struct clut_state *e, float *xyz)
{
	fill_matrix_from_controls(e);
	float (*M)[4] = e->M;

	// clip coordinates (cx,cy,cz,cw)
	float c[4] = {
		M[0][0]*xyz[0] + M[0][1]*xyz[1] + M[0][2]*xyz[2] + M[0][3],
		M[1][0]*xyz[0] + M[1][1]*xyz[1] + M[1][2]*xyz[2] + M[1][3],
		M[2][0]*xyz[0] + M[2][1]*xyz[1] + M[2][2]*xyz[2] + M[2][3],
		M[3][0]*xyz[0] + M[3][1]*xyz[1] + M[3][2]*xyz[2] + M[3][3]
	};

	// normalized device coordinates
	float xn = c[0] / c[3];
	float yn = c[1] / c[3];
	float zn = c[2] / c[3];

	// viewport coordinates
	// NOTE: this transformation is an universal constant from within
	// this program.  We never change the viewport meaning.  We just
	// change the view/projection matrix
	ij[0] = (1 + xn) * e->w / 2;
	ij[1] = (1 - yn) * e->h / 2;

	// clip focal plane
	if (zn < 0) return false;

	// clip viewport
	if (ij[0] < 0 || ij[0] >= e->w) return false;
	if (ij[1] < 0 || ij[1] >= e->w) return false;

	return true;
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

static uint8_t rgb_black[3] = {0, 0, 0};
static uint8_t rgb_white[3] = {255, 255, 255};
static uint8_t rgb_red[3] = {255, 0, 0};
static uint8_t rgb_green[3] = {0, 255, 0};
static uint8_t rgb_dgreen[3] = {0, 155, 0};
static uint8_t rgb_blue[3] = {0, 0, 255};

static void plot_pixel(int x, int y, void *e)
{
	static char c[3] = { 0, 255, 255};
	struct FTR *f = e;
	if (x == INT_MIN && y == INT_MAX)
		for (int i = 0; i < 3; i++)
			c[i] = ((char *)e)[i];
	else
		if (insideP(f->w, f->h, x, y))
			for (int i = 0; i < 3; i++)
				f->rgb[3*(f->w*y+x)+i] = c[i];
}

static void plot_segment(struct FTR *f, float a[2], float b[2], uint8_t *c)
{
	plot_pixel(INT_MIN, INT_MAX, c);
	traverse_segment(a[0], a[1], b[0], b[1], plot_pixel, f);
}

static void plot_triangle_wireframe(struct FTR *f,
		float a[2], float b[2], float c[2], uint8_t *k)
{
	plot_pixel(INT_MIN, INT_MAX, k);
	traverse_segment(a[0], a[1], b[0], b[1], plot_pixel, f);
	traverse_segment(b[0], b[1], c[0], c[1], plot_pixel, f);
	traverse_segment(c[0], c[1], a[0], a[1], plot_pixel, f);
}


// CALLBACK : expose
static void event_expose(struct FTR *f, int ev_b, int ev_m, int ev_x, int ev_y)
{
	struct clut_state *e = f->userdata;

	// light green background
	for (int i = 0; i < f->w * f->h; i++)
	{
		f->rgb[3*i+0] = 200;
		f->rgb[3*i+1] = 255;
		f->rgb[3*i+2] = 200;
	}

	// expose triangles
	//
	// 1st version: transparent wireframe
	for (int i = 0; i < e->nt; i++)
	{
		int *t = e->t[i];
		float *v[3] = { e->v[t[0]], e->v[t[1]], e->v[t[2]] };
		float V[3][2];
		for (int j = 0; j < 3; j++)
			project(V[j], e, v[j]);
		plot_triangle_wireframe(f, V[0], V[1], V[2], rgb_black);
	}

	// expoes axes
	float ax[3][2][3] = {
		{ {0,0,0}, {1,0,0} },
		{ {0,0,0}, {0,1,0} },
		{ {0,0,0}, {0,0,1} }
	};
	float A[3][2][2];
	for (int i = 0; i < 3; i++)
	for (int k = 0; k < 2; k++)
		project(A[i][k], e, ax[i][k]);
	plot_segment(f, A[0][0], A[0][1], rgb_red);
	plot_segment(f, A[1][0], A[1][1], rgb_dgreen);
	plot_segment(f, A[2][0], A[2][1], rgb_blue);


	// expose HUD
	char buf[0x200] = {0};
	snprintf(buf, 0x200,
			"a = %g\n"
			"b = %g\n"
			"c = %g\n"
			, e->a, e->b, e->c);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			0, 0+0, rgb_dgreen, rgb_black, 0, e->font, buf);


	f->changed = 1;
}


// CALLBACK: mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct clut_state *e = f->userdata;


//	// drag WINDOW DOMAIN background (realtime feedback)
//	if (e->dragging_background && m & FTR_BUTTON_LEFT)
//	{
//		int dx = x - e->drag_handle[0];
//		int dy = y - e->drag_handle[1];
//		change_view_offset(e, dx, dy);
//		e->drag_handle[0] = x;
//		e->drag_handle[1] = y;
//		f->changed = 1;
//		return;
//	}
}


// CALLBACK : resize
static void event_resize(struct FTR *f, int b, int m, int x, int y)
{
	struct clut_state *e = f->userdata;
	fprintf(stderr, "resize %d %d\n", x, y);
	e->w = f->w;
	e->h = f->h;
	f->changed = 1;
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

//static void action_offset_q0(struct jmg_state *e, float d[2])
//{
//	// get ij of current point
//	float ij[2];
//	win_from_xy(ij, e, e->x);
//
//	// shift ij by the required amount
//	ij[0] += d[0] * e->w / 2;
//	ij[1] += d[1] * e->h / 2;
//
//	// back to plane domain
//	xy_from_win(e->x, e, ij);
//}

// CALLBACK : key
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);

	//fprintf(stderr, "\tevent KEY='%c' (%d)\n", k, k);

	struct clut_state *e = f->userdata;

//	// "hidden" keys (not visible directly in the hud)
//	if (k == ',') action_screenshot(f);
//	if (k == 'm' || k == ' ') cycle_int(&e->bg_mode, 1, 4);
//	if (k == 'M' || k == '\b') cycle_int(&e->bg_mode, -1, 4);
//	if (k == 'h') scale_float(&e->nskip, 2);
//	if (k == 'H') scale_float(&e->nskip, 0.5);
//	if (k == 'd') scale_float(&e->gstep, cbrt(2));
//	if (k == 'D') scale_float(&e->gstep, 1/cbrt(2));
//	if (k == 'u') cycle_int(&e->hud, 1, 2);
//	if (e->nskip < 1) e->nskip = 1;
//	if (tolower(k)=='d')fprintf(stderr,"gstep=%g\n", e->gstep);
//
//	// same letter as they appear on the hud
	if (k == 'a') shift_float(&e->a, -5);
	if (k == 'A') shift_float(&e->a, +5);
	if (k == 'b') shift_float(&e->b, -5);
	if (k == 'B') shift_float(&e->b, +5);
	if (k == 'c') shift_float(&e->c, -5);
	if (k == 'C') shift_float(&e->c, +5);
//	if (k == 'e') shift_float(&e->E, -0.125);
//	if (k == 'E') shift_float(&e->E, +0.125);
//	if (k == 'b') shift_float(&e->bg_A, -0.125);
//	if (k == 'B') shift_float(&e->bg_A, +0.125);
//	if (k == 'j') shift_angle(&e->j0, +10);
//	if (k == 'J') shift_angle(&e->j0, -10);
//	if (k == 'n') scale_int(&e->N, 1.0/1.1);
//	if (k == 'N') scale_int(&e->N, 1.1/1.0);
//	if (k == 't') scale_float(&e->tstep, pow(2,-0.25));
//	if (k == 'T') scale_float(&e->tstep, pow(2,+0.25));
//	if (k == 's') cycle_int(&e->solver, +1, 4);
//	if (k == 'S') cycle_int(&e->solver, -1, 4);
//
//	// arrows to move the starting point
//	float d[2] = {0, 0}, inc = 0.02;
//	if (m & FTR_MASK_SHIFT  ) inc /= 10;
//	if (m & FTR_MASK_CONTROL) inc *= 10;
//	if (k == FTR_KEY_LEFT   ) d[0] -= inc;
//	if (k == FTR_KEY_RIGHT  ) d[0] += inc;
//	if (k == FTR_KEY_UP     ) d[1] -= inc;
//	if (k == FTR_KEY_DOWN   ) d[1] += inc;
//	action_offset_q0(e, d);

	f->changed = 1;
}

// CALLBACK : mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct clut_state *e = f->userdata;
	//printf("event button k=%d m=%d x=%d y=%d\n", k, m, x, y);

//	// right-click : move query point
//	if (k == FTR_BUTTON_RIGHT)
//	{
//		float ij[2] = {x, y};
//		xy_from_win(e->x, e, ij);
//	}


//	// wheel : change written parameters
//	// a, E, m, A, j0, v0, solver, N, h, Tn, Ts
//	// 0  1  2  3  4   5   6       7  8  9   10
//	// (hitboxes of font height)
//	int Y = y / e->font->height;
//	int X = x / e->font->width;
//	if (k == FTR_BUTTON_DOWN && e->hud
//			&& x < 30*e->font->width
//			&& y < 10*e->font->height
//	   )
//	{
//		if (Y == 0) shift_float(&e->a, -0.125);
//		if (Y == 1) shift_float(&e->E, -0.125);
//		if (Y == 2) shift_float(&e->bg_A, -0.125);
//		if (Y == 3) shift_angle(&e->j0, -10);
//		if (Y == 4) scale_float(&e->v0, 1/pow(2,0.25));
//		if (Y == 5) cycle_int(&e->solver, -1, 4);
//		if (Y == 6) scale_int(&e->N, 1.0/1.1);
//		if (Y == 7) scale_float(&e->tstep, 1/pow(2,0.25));
//		if (Y == 8) shift_int(&e->tissot_n, -1);
//		if (Y == 9) scale_float(&e->tissot_scale, 1/1.1);
//		f->changed = 1;
//		return;
//	}
//	if (k == FTR_BUTTON_UP && e->hud
//			&& x < 30 * e->font->width
//			&& y < 10*e->font->height
//			)
//	{
//		if (Y == 0) shift_float(&e->a, 0.125);
//		if (Y == 1) shift_float(&e->E, 0.125);
//		if (Y == 2) shift_float(&e->bg_A, 0.125);
//		if (Y == 3) shift_angle(&e->j0, 10);
//		if (Y == 4) scale_float(&e->v0, pow(2,0.25));
//		if (Y == 5) cycle_int(&e->solver, 1, 4);
//		if (Y == 6) scale_int(&e->N, 1.1);
//		if (Y == 7) scale_float(&e->tstep, pow(2,0.25));
//		if (Y == 8) shift_int(&e->tissot_n, 1);
//		if (Y == 9) scale_float(&e->tissot_scale, 1.1);
//		f->changed = 1;
//		return;
//	}


//	// begin dragging a the WINDOW BACKGROUND
//	if (k == FTR_BUTTON_LEFT)// && hit_point(e, x, y) < 0)
//	{
//		e->drag_handle[0] = x;
//		e->drag_handle[1] = y;
//		e->dragging_background = true;
//	}
//
//	// end dragging the WINDOW BACLGROUND
//	if (e->dragging_background && k == -FTR_BUTTON_LEFT)
//	{
//		int dx = x - e->drag_handle[0];
//		int dy = y - e->drag_handle[1];
//		change_view_offset(e, dx, dy);
//		e->dragging_background = false;
//	}
//
//
//	// radius in/out (if hit), zoom in/out (if no hit)
//	if (k == FTR_BUTTON_DOWN)
//		change_view_scale(e, x, y, ZOOM_FACTOR);
//	if (k == FTR_BUTTON_UP)
//		change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);

	f->changed = 1;
}

//static void print_args(int c, char **v, char *s)
//{
//	for (int i = 0; i <= c; i++)
//		printf("%s\t: ARG[%d/%d] = \"%s\"\n", s, i, c, v[i]);
//	printf("\n");
//}

#include "pickopt.c"
int main_clut(int c, char *v[])
{
	char *output_file = pick_option(&c, &v, "o", "");

	int w      = atoi(pick_variable(&c, &v, "w", "800"));
	int h      = atoi(pick_variable(&c, &v, "h", "800"));

	struct clut_state e[1];
	clut_emtpy(e);
	e->w = w;
	e->h = h;
	clut_fill_simplex(e);
	//e->a       = atof(pick_variable(&c, &v, "a", "-1"));
	//e->E       = atof(pick_variable(&c, &v, "E", "-1"));
	//e->bg_mode = atoi(pick_variable(&c, &v, "mode", "1"));
	//e->bg_A    = atof(pick_variable(&c, &v, "bg_A", "1"));
	//e->x[0]    = atof(pick_variable(&c, &v, "x", "0.7"));
	//e->x[1]    = atof(pick_variable(&c, &v, "y", "0.0"));
	//e->j0      = atof(pick_variable(&c, &v, "j", "120"));
	//e->solver  = atoi(pick_variable(&c, &v, "solver", "1"));
	//e->N       = atoi(pick_variable(&c, &v, "N", "500"));
	//e->tstep   = atof(pick_variable(&c, &v, "tstep", "0.0078125"));
	//e->nskip   = atoi(pick_variable(&c, &v, "nskip", "1"));
	//e->gstep   = atof(pick_variable(&c, &v, "gstep", "0.05"));
	//e->hud     = atoi(pick_variable(&c, &v, "hud", "1"));
	//e->tissot_n = atoi(pick_variable(&c, &v, "tissot_n", "27"));
	//e->tissot_scale = atof(pick_variable(&c, &v, "tissot_scale", "7"));

//	if (*output_file) // non-interactive mode
//	{
//		struct FTR f = {.w = w, .h = h, .userdata = e};
//		f.rgb = malloc(3*w*h);
//		event_expose(&f, 0,0,0,0);
//#ifndef __EMSCRIPTEN__
//		void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
//		iio_write_image_uint8_vec(output_file, f.rgb, w, h, 3);
//#endif
//		return 0;
//	}


	struct FTR f = ftr_new_window(e->w, e->h);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "key", event_key);
	ftr_loop_run(&f);
	ftr_close(&f);
	return 0;
}

#ifdef __EMSCRIPTEN__
int main(void)
{
	int c = 1;
	char *v[2] = {"clut", NULL};
#else
int main(int c, char *v[])
{
#endif
	return main_clut(c, v);
}
