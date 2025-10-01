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
#define MAX_KAPPAS 20

// radius of the points
#define POINT_RADIUS 5.3
//#define POINT_RADIUS 1

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
	int mode;             // LEFT,RIGHT,MIN,MAX,ALL [0..4]
	int action;           // EQUAL,EXPONENTIAL,RANDOM,JITTER

	// display theme
	uint8_t rgb_bg[3];
	uint8_t rgb_fg[3];
	uint8_t rgb_axes[3];
	uint8_t rgb_grid[3];
	uint8_t rgb_curv[3];
	uint8_t rgb_nan[3];
	struct bitmap_font fonts[5]; // from small to large
	struct bitmap_font *font;

};


// reset and center the state
static void init_state(struct viewer_state *e)
{
	e->n = 2;
	e->k[0] = 0.3;
	e->k[1] = 0.7;
	e->mode = 4;
	//e->action = 1;
	//e->m = 1;
	//e->n = 1;
	//e->q = 3;


	e->w = 800;
	e->h = 800;
	e->x0 = -0.2;
	e->xf = 1.2;
	e->y0 = -0.2;
	e->yf = 1.2;

	e->rgb_bg[0] = 0;
	e->rgb_bg[1] = 0;
	e->rgb_bg[2] = 0;

	e->rgb_fg[0] = 100;
	e->rgb_fg[1] = 200;
	e->rgb_fg[2] = 150;

	e->rgb_axes[0] = 100;
	e->rgb_axes[1] = 100;
	e->rgb_axes[2] = 100;

	e->rgb_grid[0] = 40;
	e->rgb_grid[1] = 40;
	e->rgb_grid[2] = 40;

	e->rgb_curv[0] = 255;
	e->rgb_curv[1] = 100;
	e->rgb_curv[2] = 0;

	e->font = e->fonts + 3;

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
	{0, 255, 0},      // 10  green (go for a "oscilloscope" aesthetic)
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

static void get_win_from_xy(struct viewer_state *e, float ij[2], float xy[2])
{
	ij[0] =            e->w * (xy[0] - e->x0) / (e->xf - e->x0);
	ij[1] = e->h - 1 - e->h * (xy[1] - e->y0) / (e->yf - e->y0);
}

static int compare_points_lexicographically(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	int p = (a[0] > b[0]) - (a[0] < b[0]);
	if (p)
		return p;
	else
		return ((a[1] > b[1]) - (a[1] < b[1]));
}

// NOTE: T must hold 2^(n+1) floats
static void sorted_kappa_sums(float *T, struct viewer_state *e)
{
	int N = 1 << e->n;
	for (int i = 0; i < N; i++)
	{
		float k = 0; // accumulator for this kappa sum
		for (int j = 0; j < e->n; j++)
			if (i & (1 << j))
				k += e->k[j];
		T[2*i+0] = k; // kappa sum
		T[2*i+1] = i; // kappas to be summed (in the bits of i)
	}
	qsort(T, N, 2*sizeof*T, compare_points_lexicographically);
}

static float kappa_sum(struct viewer_state *e)
{
	float k = 0;
	for (int i = 0; i < e->n; i++)
		k += e->k[i];
	return k;
}

static void assert_kappa_sum(struct viewer_state *e)
{
	float k = kappa_sum(e);
	if (fabs(k - 1) > 1e-5)
	{
		fprintf(stderr, "ksum = %g (diff=%g)\n", k, k-1);
		assert(false);
	}
}

static void set_equal_kappas(struct viewer_state *e)
{
	for (int i = 0; i < e->n; i++)
		e->k[i] = 1.0 / e->n;
	assert_kappa_sum(e);
}

static void set_exponential_kappas(struct viewer_state *e)
{
	for (int i = 0; i < e->n; i++)
		e->k[i] = (1 << i) / (-1.0 + (1 << e->n));
	assert_kappa_sum(e);
}

static void set_random_kappas(struct viewer_state *e)
{
	for (int i = 0; i < e->n; i++) e->k[i] = random_uniform();
	float k = kappa_sum(e);
	for (int i = 0; i < e->n; i++) e->k[i] /= k;
	assert_kappa_sum(e);
}

static void jitter_kappas(struct viewer_state *e)
{
	for (int i = 0; i < e->n; i++)
		e->k[i] += 0.005 * random_normal();
	for (int i = 0; i < e->n; i++)
		e->k[i] = fmax(0, fmin(1, e->k[i]));
	float k = kappa_sum(e);
	for (int i = 0; i < e->n; i++) e->k[i] /= k;
	assert_kappa_sum(e);
}

static void set_new_n(struct viewer_state *e, int n)
{
	if (n < 1 || n > MAX_KAPPAS)
		return;
	e->n = n;
	set_exponential_kappas(e);
}

static void change_one_kappa(struct viewer_state *e, int i, float d)
{
	assert(i >= 0 && i < e->n);
	switch (e->mode) {
	case 0: // LEFT
		break;
	case 1: // RIGHT
		break;
	case 2: // MIN
		break;
	case 3: // MAX
		break;
	case 4: // ALL
		for (int j = 0; j < e->n; j++)
			if (j == i)
				e->k[j] += d;
			else
				e->k[j] -= d/(e->n - 1);
		break;
	}
}




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

static void plot_pixel(int x, int y, void *e)
{
	static char c[3] = { 0, 255, 255};
	struct FTR *f = e;
	if (x == -1 && y == -1)
		for (int i = 0; i < 3; i++)
			c[i] = ((char *)e)[i];
	else
		if (insideP(f, x, y))
			for (int i = 0; i < 3; i++)
				f->rgb[3*(f->w*y+x)+i] = c[i];
}

static void plot_segment(struct FTR *f,
		float x0, float y0, float xf, float yf, uint8_t c[3])
{
	plot_pixel(-1, -1, c);
	traverse_segment(x0, y0, xf, yf, plot_pixel, f);
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


static void plot_line_in_plane(struct FTR *f, float *P, float *Q, uint8_t c[3])
{
	struct viewer_state *e = f->userdata;
	float p[2], q[2];
	get_win_from_xy(e, p, P);
	get_win_from_xy(e, q, Q);
	plot_segment(f, p[0], p[1], q[0], q[1], c);
}

static int count_bits(unsigned int n)
{
	int o = n;
	int c = 0;
	for (; n; c++)
		n &= n - 1;  // clear the least significant bit set
	return c;
}

static long binomial(int n, int k)
{
	static int T[MAX_KAPPAS+1][MAX_KAPPAS+1] = {0};
	int N = MAX_KAPPAS;
	if (!T[0][0])
	{
		for (int i = 0; i <= N; i++)
			T[i][0] = T[i][i] = 1;
		for (int i = 2; i <= N; i++)
		for (int j = 1; j < i; j++)
			T[i][j] = T[i-1][j-1] + T[i-1][j];
		//for (int i = 0; i <= N; i++)
		//for (int j = 0; j <= i; j++)
		//	fprintf(stderr, "T[%d][%d] = %d\n", i, j, T[i][j]);
	}
	assert(n > 0 && n <= N);
	assert(k >= 0 && k <= N);
	return T[n][k];
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
	for (int x = floor(e->x0); x < ceil(e->xf); x++)
	for (int j = 0; j < e->h; j++)
	{
		int i = e->w * (x - e->x0) / (e->xf - e->x0);
		uint8_t *c = x ? e->rgb_grid : e->rgb_axes;
		if (insideP(f, i, j))
		for (int k = 0; k < 3; k++)
			f->rgb[(f->w*j+i)*3+k] = c[k];
	}

	// plot grid: horizontal lines
	for (int y = floor(e->y0); y < ceil(e->yf); y++)
	for (int i = 0; i < e->w; i++)
	{
		int j = e->h - 1 - e->h * (y - e->y0) / (e->yf - e->y0);
		uint8_t *c = y ? e->rgb_grid : e->rgb_axes;
		if (insideP(f, i, j))
		for (int k = 0; k < 3; k++)
			f->rgb[(f->w*j+i)*3+k] = c[k];
	}

	float ab[4][2] = { {0.5, 0}, {0.5, 1}, {0, 0.5}, {1, 0.5} };
	float AB[4][2];
	for (int i = 0; i < 4; i++)
		get_win_from_xy(e, AB[i], ab[i]);
	plot_segment(f, AB[0][0], AB[0][1], AB[1][0], AB[1][1], e->rgb_grid);
	plot_segment(f, AB[2][0], AB[2][1], AB[3][0], AB[3][1], e->rgb_grid);


	//// plot kappas
	//for (int i = 0; i < e->n; i++)
	//{
	//	float xy[2] = {e->k[i], 0}, ij[2];
	//	get_win_from_xy(e, ij, xy);
	//	splat_disk(f->rgb, f->w, f->h, ij, POINT_RADIUS, e->rgb_curv);
	//}

	// plot kappa sums
	int N = 1 << e->n; // number of kappa sums = 2^n
	for (int i = 0; i < N; i++)
	{
		float k = 0; // accumulator for this kappa sum
		for (int j = 0; j < e->n; j++)
			if (i & (1 << j))
				k += e->k[j];
		float r = (i&&!(i&(i-1))) ? 7.3 : 2.7; // pure kappa big dot
		float xy[2] = {k, 0}, ij[2];
		get_win_from_xy(e, ij, xy);
		splat_disk(f->rgb, f->w, f->h, ij, r, e->rgb_curv);
	}

	// diagonal line (accumulated uniform distribution)
	plot_line_in_plane(f, (float[]){0,0}, (float[]){1,1}, palette[3]);


	// bound lines
	float M = -INFINITY;
	for (int i = 0; i < e->n; i++)
		M = fmax(M, e->k[i]);
	fprintf(stderr, "M = %g\n", M);
	plot_line_in_plane(f, (float[]){M,0}, (float[]){1,1-M}, palette[5]);
	plot_line_in_plane(f, (float[]){0,M}, (float[]){1-M,1}, palette[5]);

	// build table of kappa sums
	float T[2*N];
	sorted_kappa_sums(T, e);

	// plot accumulated distribution
	float P[2] = {e->x0, 0}, *p = P;
	float Q[2] = {0    , 0}, *q = Q;
	plot_line_in_plane(f, p, q, palette[14]);
	for (int i = 1; i < N; i++)
	{
		p[0] = q[0];
		//p[1] = q[1] + 1.0 / (1 << e->n);
		int k = count_bits(i);
		p[1] = q[1] + 1.0 /((e->n + 1.0) * binomial(e->n, k));
		// TODO: fix the formula on previous line
		q[0] = T[2*i+0];
		q[1] = p[1];
		plot_line_in_plane(f, p, q, palette[14]);
	}
	p[0] = 1   ;  p[1] = 1;
	q[0] = e->xf; q[1] = 1;
	plot_line_in_plane(f, p, q, palette[14]);


	//for (int i = 0; i < e->w; i++)
	//{
	//	float x = e->x0 + i * (e->xf - e->x0) / e->w;
	//	float y = F[e->m](e->n, e->q, x);
	//	float j = e->h - 1 - e->h * (y - e->y0) / (e->yf - e->y0);
	//	float P[2] = {i, j};
	//	//if (0 == i%10)
	//	splat_disk(f->rgb, f->w, f->h, P, POINT_RADIUS, e->rgb_curv);
	//}


	// hud
	char buf[0x400];
	int b = 0;
	b += snprintf(buf+b, 0x400-b, "n = %d\nk = ", e->n);
	for (int i = 0; i < e->n; i++)
		b += snprintf(buf+b, 0x400-b, "%g\t", e->k[i]);
	b += snprintf(buf+b, 0x400-b, "\nmode = %d LEFT, RIGHT, MIN, MAX, ALL\n", e->mode);
	b += snprintf(buf+b, 0x400-b, "action = UNI, EXP, RND, JIT");
	put_string_in_rgb_image(f->rgb, f->w, f->h,
			10, 0, e->rgb_fg, e->rgb_bg, 0, e->font, buf);
}



// SECTION 8. User-Interface Actions and Events                             {{{1

// action: change n
static void change_n(struct viewer_state *e, int d)
{
	int n = e->n + d;
	if (n < 1) n = 1;
	if (n > MAX_KAPPAS) n = MAX_KAPPAS;
	set_new_n(e, n);
}

// action: viewport translation
static void change_view_offset(struct viewer_state *e, float dx, float dy)
{
	e->offset[0] += dx;
	e->offset[1] += dy;
}

//// action: viewport zoom
//static void change_view_scale(struct viewer_state *e, int x, int y, float fac)
//{
//	float center[2], X[2] = {x, y};
//	map_window_to_view(e, center, X);
//	e->scale *= fac;
//	for (int p = 0; p < 2; p++)
//		e->offset[p] = -center[p]*e->scale + X[p];
//	fprintf(stderr, "zoom changed %g\n", e->scale);
//}

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
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);

	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'n') change_n(e, 1);
	if (k == 'N') change_n(e, -1);
	if (k == 'j') jitter_kappas(e);
	if (k == 'r') set_random_kappas(e);
	if (k == 'e') set_equal_kappas(e);
	if (k == 'x') set_exponential_kappas(e);
	if (k == ',') action_screenshot(f);

	//e->dragging_window_point = e->dragging_image_point = false;
	f->changed = 1;
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	f->changed = 1;
}

//static void toggle_mode(struct viewer_state *e, int s)
//{
//	e->m += s;
//	if (e->m < 0) e->m = 0;
//	if (e->m >= 6) e->m = 5;
//	fprintf(stderr, "e->m = %d\n", e->m);
//}
//static void shift_param_n(struct viewer_state *e, float s) { e->n += s; }
//static void shift_param_q(struct viewer_state *e, float s) { e->q += s; }

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

//	// m, n, q hitboxes of font height
//	// 0  1  2
//	int Y = y / e->font->height;
//	if (k == FTR_BUTTON_DOWN)
//	{
//		if (Y == 0) toggle_mode(e, -1);
//		if (Y == 1) shift_param_n(e, -1);
//		if (Y == 2) shift_param_q(e, -1);
//	}
//	if (k == FTR_BUTTON_UP)
//	{
//		if (Y == 0) toggle_mode(e, +1);
//		if (Y == 1) shift_param_n(e, +1);
//		if (Y == 2) shift_param_q(e, +1);
//	}

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

	//// drag WINDOW DOMAIN background (realtime feedback)
	//if (e->dragging_background && m & FTR_BUTTON_LEFT)
	//{
	//	int dx = x - e->drag_handle[0];
	//	int dy = y - e->drag_handle[1];
	//	change_view_offset(e, dx, dy);
	//	e->drag_handle[0] = x;
	//	e->drag_handle[1] = y;
	//	f->changed = 1;
	//	return;
	//}

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
	e->fonts[3] = reformat_font(*xfont_9x18B, UNPACKED);//the only one used
	e->fonts[4] = reformat_font(*xfont_10x20, UNPACKED);

	// open the window
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
