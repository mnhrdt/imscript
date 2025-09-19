// interactive visualization of kappa sums (for the KRT statistics)

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

// maximum allowed number of kappa points
#define MAX_KAPPAS

// radius of the points
//#define POINT_RADIUS 2.3
#define POINT_RADIUS 1

// zoom factor for zoom-in and zoom-out
#define ZOOM_FACTOR 1.43

// radius scaling factor for the inversion circle
#define RADIUS_FACTOR 1.13

// data structure to store the state of the viewer
struct viewer_state {

	// view port geometry
	int w, h;              // view field dimensions in pixels
	float x0, xf, y0, yf;  // view rectangle on R^2
	float offset[2], scale;  // view rectangle on R^2

	// visualization parameters
	int n;                // number of kappas (1,2,...,MAX_KAPPAS)
	float k[MAX_KAPPAS];  // values of each kappa
	int mode;             // -1=put_left 0=put_all 1=put_right
	int preset;           // 0=uniform 1=random

	// display theme
	uint8_t rgb_bg[3];
	uint8_t rgb_fg[3];
	uint8_t rgb_axes[3];
	uint8_t rgb_grid[3];
	uint8_t rgb_curv[3];
	uint8_t rgb_nan[3];
	struct bitmap_font fonts[5]; // from small to large
	struct bitmap_font *font;

	// drag state
	bool dragging_window_point;
	bool dragging_image_point;
	bool dragged_point;
	bool dragging_background;
	int drag_handle[2];
};


// reset and center the state
static void init_state(struct viewer_state *e)
{
	//e->m = 1;
	//e->n = 1;
	//e->q = 3;


	e->w = 800;
	e->h = 600;
	e->x0 = -1;
	e->xf = 7;
	e->y0 = -3;
	e->yf = 3;

	e->n = 3;
	e->mode = 1;
	e->preset = 0;

	e->rgb_bg[0] = 0;
	e->rgb_bg[1] = 0;
	e->rgb_bg[2] = 0;

	e->rgb_fg[0] = 100;
	e->rgb_fg[1] = 200;
	e->rgb_fg[2] = 150;

	e->rgb_axes[0] = 200;
	e->rgb_axes[1] = 200;
	e->rgb_axes[2] = 200;

	e->rgb_grid[0] = 40;
	e->rgb_grid[1] = 40;
	e->rgb_grid[2] = 40;

	e->rgb_curv[0] = 255;
	e->rgb_curv[1] = 100;
	e->rgb_curv[2] = 0;

	e->font = e->fonts + 3;

	// drag state
	e->dragging_window_point = false;
	e->dragging_image_point = false;
	e->dragged_point = -1;
	e->dragging_background = false;
}


// xterm-like 16 color palette
static uint8_t palette[16][3] = {
	{0, 0, 0},        //  0  black
	{128, 0, 0},      //  1  dark red
	{0, 128, 0},      //  2  dark green
	{128, 128, 0},    //  3  dark yellow
	{0, 0, 128},      //  4  dark blue
	{128, 0, 128},    //  5  dark magenta
	{0, 128, 128},    //  6  dark cyan
	{192, 192, 192},  //  7  dark white
	{128, 128, 128},  //  8  gray
	{255, 0, 0},      //  9  red
	{0, 255, 0},      // 10  green
	{255, 255, 0},    // 11  yellow
	{0, 0, 255},      // 12  blue
	{255, 0, 255},    // 13  magenta
	{0, 255, 255},    // 14  cyan
	{255, 255, 255},  // 15  white
};


// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 3. algorithms                                                    {{{1






// SECTION 7. Drawing                                                       {{{1

// Subsection 7.1. Drawing segments                                         {{{2

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

// draw a segment between two points
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

// auxiliary function for drawing a red pixel
static void plot_pixel_red(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}

// auxiliary function for drawing a blue pixel
static void plot_pixel_blue(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 255;
	}
}

// auxiliary function for drawing a cyan pixel
static void plot_pixel_cyan(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 255;
		f->rgb[3*idx+2] = 255;
	}
}

// auxiliary function for drawing a green pixel
static void plot_pixel_green(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 128;
		f->rgb[3*idx+2] = 0;
	}
}

// auxiliary function for drawing a gray pixel
static void plot_pixel_gray(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 120;
		f->rgb[3*idx+1] = 120;
		f->rgb[3*idx+2] = 120;
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

// function to draw a blue segment
static void plot_segment_blue(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_blue, f);
}

// function to draw a cyan segment
static void plot_segment_cyan(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_cyan, f);
}

// function to draw a gray segment
static void plot_segment_gray(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_gray, f);
}

// function to draw a green segment
static void plot_segment_green(struct FTR *f,
		float x0, float y0, float xf, float yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_green, f);
}

static void plot_circle_green(struct FTR *f,
		float x, float y, float r)
{
	traverse_circle(x, y, r, plot_pixel_green, f);
}



// Subsection 7.2. Drawing user-interface elements                          {{{2

static void splat_disk(uint8_t *rgb, int w, int h, float p[2], float r,
		uint8_t color[3])
{
	for (int j = -r-1 ; j <= r+1; j++)
	for (int i = -r-1 ; i <= r+1; i++)
	if (hypot(i, j) < r)
	{
		int ii = p[0] + i;
		int jj = p[1] + j;
		if (ii>=0 && jj>=0 && ii<w && jj<h)
		//{
		//	float a = pow(hypot(i, j)/r, 4);
			for (int k = 0; k < 3; k++)
				//rgb[3*(w*jj+ii)+k] = a*254+(1-a)*color[k];
		//		rgb[3*(w*jj+ii)+k] = a*0+(1-a)*color[k];
				rgb[3*(w*jj+ii)+k] = color[k];
		//}
	}
}





// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// clear canvas
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = e->rgb_bg[i%3];

	// plot grid: vertical lines
	for (int x = ceil(e->x0); x < floor(e->xf); x++)
	for (int j = 0; j < e->h; j++)
	{
		int i = e->w * (x - e->x0) / (e->xf - e->x0);
		uint8_t *c = x ? e->rgb_grid : e->rgb_axes;
		if (insideP(f, i, j))
		for (int k = 0; k < 3; k++)
			f->rgb[(f->w*j+i)*3+k] = c[k];
	}

	// plot grid: horizontal lines
	for (int y = ceil(e->y0); y < floor(e->yf); y++)
	for (int i = 0; i < e->w; i++)
	{
		int j = e->h - 1 - e->h * (y - e->y0) / (e->yf - e->y0);
		uint8_t *c = y ? e->rgb_grid : e->rgb_axes;
		if (insideP(f, i, j))
		for (int k = 0; k < 3; k++)
			f->rgb[(f->w*j+i)*3+k] = c[k];
	}

	for (int i = 0; i < e->w; i++)
	{
		float x = e->x0 + i * (e->xf - e->x0) / e->w;
		float y = F[e->m](e->n, e->q, x);
		float j = e->h - 1 - e->h * (y - e->y0) / (e->yf - e->y0);
		float P[2] = {i, j};
		//if (0 == i%10)
		splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, e->rgb_curv);
	}

	// hud
	char buf[0x400];
	char *M[6] = {"ce", "se", "Mc1", "Ms1", "Mc2", "Ms2"};
	snprintf(buf, 0x400, "m = %s\nn = %g\nq = %g\n", M[e->m], e->n, e->q);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			10, 0, e->rgb_fg, e->rgb_bg, 0, e->font, buf);
}



// SECTION 8. User-Interface Actions and Events                             {{{1


// action: viewport translation
static void change_view_offset(struct viewer_state *e, float dx, float dy)
{
	e->offset[0] += dx;
	e->offset[1] += dy;
}

// action: viewport zoom
static void change_view_scale(struct viewer_state *e, int x, int y, float fac)
{
	float center[2], X[2] = {x, y};
	map_window_to_view(e, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "zoom changed %g\n", e->scale);
}

static void action_screenshot(struct FTR *f)
{
	static int c = 0;
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "screenshot_cloudette_%d.png", c);
	void iio_write_image_uint8_vec(char*,uint8_t*,int,int,int);
	iio_write_image_uint8_vec(n, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}


// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'j') change_view_offset(e, 0, -10);
	if (k == 'k') change_view_offset(e, 0, 10);
	if (k == 'h') change_view_offset(e, 10, 0);
	if (k == 'l') change_view_offset(e, -10, 0);
	if (k == FTR_KEY_DOWN ) change_view_offset(e, 0, -100);
	if (k == FTR_KEY_UP   ) change_view_offset(e, 0, 100);
	if (k == FTR_KEY_RIGHT) change_view_offset(e, -100, 0);
	if (k == FTR_KEY_LEFT)  change_view_offset(e, 100, 0);
	if (k == '+') change_view_scale(e, f->w/2, f->h/2, ZOOM_FACTOR);
	if (k == '-') change_view_scale(e, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);
	//if (k == 'a') e->restrict_to_affine = !e->restrict_to_affine;
	if (k == ',') action_screenshot(f);

	e->dragging_window_point = e->dragging_image_point = false;
	f->changed = 1;
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	f->changed = 1;
}

static void toggle_mode(struct viewer_state *e, int s)
{
	e->m += s;
	if (e->m < 0) e->m = 0;
	if (e->m >= 6) e->m = 5;
	fprintf(stderr, "e->m = %d\n", e->m);
}
static void shift_param_n(struct viewer_state *e, float s) { e->n += s; }
static void shift_param_q(struct viewer_state *e, float s) { e->q += s; }

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// m, n, q hitboxes of font height
	// 0  1  2
	int Y = y / e->font->height;
	if (k == FTR_BUTTON_DOWN)
	{
		if (Y == 0) toggle_mode(e, -1);
		if (Y == 1) shift_param_n(e, -1);
		if (Y == 2) shift_param_q(e, -1);
	}
	if (k == FTR_BUTTON_UP)
	{
		if (Y == 0) toggle_mode(e, +1);
		if (Y == 1) shift_param_n(e, +1);
		if (Y == 2) shift_param_q(e, +1);
	}

	f->changed = 1;
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;


	//// drag WINDOW DOMAIN control point (realtime feedback)
	//if (e->dragging_window_point && m & FTR_BUTTON_LEFT)
	//{
	//	drag_point_in_window_domain(e, x, y);
	//	f->changed = 1;
	//}

	// drag WINDOW DOMAIN background (realtime feedback)
	if (e->dragging_background && m & FTR_BUTTON_LEFT)
	{
		int dx = x - e->drag_handle[0];
		int dy = y - e->drag_handle[1];
		change_view_offset(e, dx, dy);
		e->drag_handle[0] = x;
		e->drag_handle[1] = y;
		f->changed = 1;
		return;
	}

//	int p = hit_point(e, x, y);
//	if (p >= 0)
//	{
//		//fprintf(stderr, "hit point %d (%d %d)\n", p, x, y);
//		e->p = p;
//		e->px = x;
//		e->py = y;
//		f->changed = 1;
//	}
//	if (p < 0 && e->p >= 0)
//	{
//		e->p = -1;
//		f->changed = 1;
//	}

}

// expose handler
static void event_expose(struct FTR *f, int b, int m, int x, int y)
{
	if (f->changed)
		paint_state(f);
}




// SECTION 10. Main Program                                                 {{{1

int main_kappaview(int argc, char *argv[])
{
	if (argc != 1)
		return fprintf(stderr, "usage:\n\t%s\n", *argv);


	// initialize state with a default view
	struct viewer_state e[1];
	init_state(e);

	// init fonts
	e->fonts[0] = reformat_font(*xfont_4x6, UNPACKED);
	e->fonts[1] = reformat_font(*xfont_6x12, UNPACKED);
	e->fonts[2] = reformat_font(*xfont_7x13, UNPACKED);
	e->fonts[3] = reformat_font(*xfont_9x18B, UNPACKED);//used
	e->fonts[4] = reformat_font(*xfont_10x20, UNPACKED);

	// open the window
	//struct FTR f = ftr_new_window(800,600);
	struct FTR f = ftr_new_window(e->w, e->h);
	f.userdata = e;
	f.changed = 1;

	// set event handlers
	ftr_set_handler(&f, "expose", event_expose);
	ftr_set_handler(&f, "resize", event_resize);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "key", event_key);

	return ftr_loop_run(&f);
}

int main(int c, char *v[]) { return main_kappaview(c, v); }

// vim:set foldmethod=marker:
