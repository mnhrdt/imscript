// cc -O2 epiview.c iio.o -o epiview -lX11 -ltiff -lpng
//
// A program for visualizing fundamental matrices
// (based in "dosdo", which was based in "vnav", which was based in "fpan")
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>

#include <stdint.h>
#include "iio.h"

#include "ftr.h"

#define WHEEL_FACTOR 1.4

#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfont_9x15.c"

// VISUALIZATION:
// * A single window split in half.
// * Left side: left image
// * Right side: right image
// * Both sides admit scroll (drag), zoom (wheel) and contrast changes
// * Hovering on the left shows the corresponding epipolar line on the right
// * And vice versa


// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. input image pair and fundamental matrix
	int w[2], h[2], pd[2];
	float *x[2];
	double fm[9], fmtr[9];

	// 2. view port parameters
	int octave[2];
	double zoom[2], offset[2][2];
	double aaa[2][3], bbb[2][3];

	// 3. window parameters
	int win_half;

	// 4. local state
	int show_lines;
	int scroll_domain;
	double dip_abc[3], dip_pq[3];
	int dip_idx, dip_xy[2];
	int lock_transform;

	// 5. visualization details
	int head_up_display;
	struct bitmap_font font;

	// TODO: allow to click on a point to fix the visualization of its line
};

// return the index of the subwindow
// transform the absolute coordinates to relative
static int subwindow(struct pan_state *e, int *i, int *j)
{
	if (*i < e->win_half)
		return 0;
	else {
		*i -= e->win_half;
		return 1;
	}
}

// change of coordinates: from window "int" pixels to image "double" point
// return the index of the image
static int window_to_image(double p[2], struct pan_state *e, int i, int j)
{
	int idx = subwindow(e, &i, &j);
	p[0] = e->offset[idx][0] + i / e->zoom[idx];
	p[1] = e->offset[idx][1] + j / e->zoom[idx];
	return idx;
}

static void image_to_window(int ij[2], struct pan_state *e, int idx,
		double x, double y)
{
	ij[0] = lrint(e->zoom[idx] * (x - e->offset[idx][0]));
	ij[1] = lrint(e->zoom[idx] * (y - e->offset[idx][1]));
	if (idx == 1)
		ij[0] += e->win_half;
}

static void matrix_transpose(double At[9], double A[9])
{
	for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
		At[i*3+j] = A[j*3+i];
}

static void matrix_times_vector(double Ax[3], double A[9], double x[3])
{
	Ax[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
	Ax[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
	Ax[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// cut a line with a segment (returns true if they cut)
static bool cut_line_with_segment(double out[2], double line[3],
		double p[2], double q[2])
{
	// points in "oriented" projective coordinates
	double pp[3] = {p[0], p[1], 1};
	double qq[3] = {q[0], q[1], 1};

	// sign of each point (says on which side of the line each point is)
	double sp = scalar_product(pp, line);
	double sq = scalar_product(qq, line);

	// if signs are different, the line crosses the segment
	if (sp * sq < 0) {
		// line trough points p and q
		double pq[3]; vector_product(pq, pp, qq);

		// intersection of "line" and "pq"
		double ii[3]; vector_product(ii, pq, line);

		// recover affine coordinates
		out[0] = ii[0] / ii[2];
		out[1] = ii[1] / ii[2];
		return true;
	}
	return false;
}

static bool cut_line_with_rectangle(double out_a[2], double out_b[2],
		double line[3], double rec_from[2], double rec_to[4])
{
	// four vertices of the rectangle
	double v[4][2] = {
		{ rec_from[0], rec_from[1] },
		{ rec_to[0]  , rec_from[1] },
		{ rec_to[0]  , rec_to[1]   },
		{ rec_from[0], rec_to[1]   }
	};

	// intersections with each of the edges
	bool xP[4]; // whether it intersects
	double x[4][2]; // where it intersects
	for (int i = 0; i < 4; i++)
		xP[i] = cut_line_with_segment(x[i], line, v[i], v[ (i+1)%4 ] );

	// write output
	int n_intersections = xP[0] + xP[1] + xP[2] + xP[3];
	if (n_intersections == 2) { // generic case: 2 intersections
		int cx = 0;
		for (int i = 0; i < 4; i++)
			if (xP[i])
			{
				double *out = cx ? out_b : out_a;
				out[0] = x[i][0];
				out[1] = x[i][1];
				cx += 1;
			}
		return true;
	}
	return false;
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

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (l >= pd) l = pd - 1;
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[pd*(j*w+i)+l];
}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	int idx = e->scroll_domain;

	e->offset[idx][0] -= dx/e->zoom[idx];
	e->offset[idx][1] -= dy/e->zoom[idx];

	if (e->lock_transform) {
		idx = !idx;
		e->offset[idx][0] -= dx/e->zoom[idx];
		e->offset[idx][1] -= dy/e->zoom[idx];
	}

	f->changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->octave[0] = e->octave[1] = 0;
	e->zoom[0] = e->zoom[1] = 1;
	e->offset[0][0] =e->offset[0][1] =e->offset[1][0] =e->offset[1][1] = 0;

	for (int i = 0; i < 2; i++)
	for (int l = 0; l < 3; l++)
	{
		e->aaa[i][l] = 1;
		e->bbb[i][l] = 0;
	}

	f->changed = 1;
}

static void action_change_zoom_to_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	int idx = subwindow(e, &x, &y);

	if (F == 1) e->octave[idx] = 0;
	e->zoom[idx] = 1/F;
	e->offset[idx][0] = p[0] - x/e->zoom[idx];
	e->offset[idx][1] = p[1] - y/e->zoom[idx];

	if (e->lock_transform) {
		idx = !idx;
		x += (idx ? -0.5 : 0.5) * e->win_half;
		if (F == 1) e->octave[idx] = 0;
		e->zoom[idx] = 1/F;
		e->offset[idx][0] = p[0] - x/e->zoom[idx];
		e->offset[idx][1] = p[1] - y/e->zoom[idx];
	}

	f->changed = 1;
}

static void action_increase_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	int xy[2] = {x, y};
	int idx = subwindow(e, xy, xy+1);

	if (e->octave[idx] < 10) {
		e->octave[idx] += 1;
		double fac = 1 << e->octave[idx];
		if (e->octave[idx] < 0) fac = 1.0/(1<<-e->octave[idx]);
		action_change_zoom_to_factor(f, x, y, fac);
	}
}

static void action_decrease_octave(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	int xy[2] = {x, y};
	int idx = subwindow(e, xy, xy+1);

	if (e->octave[idx] > 0) {
		e->octave[idx] -= 1;
		double fac = 1 << e->octave[idx];
		action_change_zoom_to_factor(f, x, y, fac);
	}
	else if (e->octave[idx] <= 0) {
		e->octave[idx] -= 1;
		double fac = 1.0/(1 << -e->octave[idx]);
		action_change_zoom_to_factor(f, x, y, fac);
	}
}

static void action_qauto(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	//float m = INFINITY, M = -m;
	float m = 0, M = 255;

	for (int i = 0; i < 2; i++)
	for (int l = 0; l < 3; l++)
	{
		e->aaa[i][l] = 255 / ( M - m );
		e->bbb[i][l] = 255 * m / ( m - M );
	}

	f->changed = 1;
}

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2], c[3];
	int i = window_to_image(p, e, x, y);
	c[0] = getsample_0(e->x[i], e->w[i], e->h[i], e->pd[i], p[0], p[1], 0);
	c[1] = getsample_0(e->x[i], e->w[i], e->h[i], e->pd[i], p[0], p[1], 1);
	c[2] = getsample_0(e->x[i], e->w[i], e->h[i], e->pd[i], p[0], p[1], 2);

	for (int l = 0; l < 3; l++)
		e->bbb[i][l] = 127.5 - e->aaa[i][l] * c[l];

	f->changed = 1;
}

static void action_contrast_span(struct FTR *f, int idx, float factor)
{
	struct pan_state *e = f->userdata;

	for (int l = 0; l < 3; l++)
	{
		float c = (127.5 - e->bbb[idx][l])/ e->aaa[idx][l];
		e->aaa[idx][l] *= factor;
		e->bbb[idx][l] = 127.5 - e->aaa[idx][l] * c;
	}

	f->changed = 1;
}

static void action_save_shot(struct FTR *f)
{
	static int shot_counter = 1;
	char fname[FILENAME_MAX];
	snprintf(fname, FILENAME_MAX, "/tmp/epiview_shot_%d.png", shot_counter);
	iio_write_image_uint8_vec(fname, f->rgb, f->w, f->h, 3);
	fprintf(stderr, "saved shot \"%s\"\n", fname);
	shot_counter += 1;
}


static void action_toggle_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->head_up_display = !e->head_up_display;
	f->changed = 1;
}

static void action_toggle_lock(struct FTR *f)
{
	struct pan_state *e = f->userdata;
	e->lock_transform = !e->lock_transform;
}

static unsigned char float_to_byte(float x)

{
	if (isnan(x)) return 0;
	if (x < 0) return 0;
	if (x > 255) return 255;
	return x;
}

static void dump_hud(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	uint8_t fg[3] = {0, 255, 0};
	uint8_t bg[3] = {0, 0, 0};
	char buf[0x100];

	snprintf(buf, 0x100, "abc");
	put_string_in_rgb_image(f->rgb,f->w,f->h,e->win_half,0,fg,bg,0,&e->font,buf);
}

static void fplot_red(int i, int j, void *ee)
{
	struct FTR *f = ee;
	if (i > 0 && j > 0 && i < f->w && j < f->h)
	{
		f->rgb[3*(j*f->w+i)+0] = 255;
		f->rgb[3*(j*f->w+i)+1] = 0;
		f->rgb[3*(j*f->w+i)+2] = 0;
	}
}

// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	for (int i = 0; i < f->w * f->h * 3; i++) f->rgb[i] = 0;

	// render both images
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		int idx = window_to_image(p, e, i, j);
		unsigned char *dest = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < 3; l++)
		{
			float v = getsample_0(e->x[idx],
					e->w[idx], e->h[idx], e->pd[idx],
					p[0], p[1], l);
			dest[l] = float_to_byte(e->aaa[idx][l] * v + e->bbb[idx][l]);
		}
	}

	// render epipolars
	if (e->show_lines && isfinite(e->dip_abc[0]))
	{
		int idx = e->dip_idx;
		// idx = 0 : dip left, line right
		// idx = 1 : dip right, line left
		double from[2], to[2], aa[2], bb[2];
		window_to_image(from, e, idx ? 0 : e->win_half, 0);
		window_to_image(to, e, idx ? e->win_half - 1 : f->w, f->h - 1);
		if (cut_line_with_rectangle(aa, bb, e->dip_abc, from, to)) {
			int iaa[2], ibb[2];
			image_to_window(iaa, e, !idx, aa[0], aa[1]);
			image_to_window(ibb, e, !idx, bb[0], bb[1]);
			traverse_segment(iaa[0], iaa[1], ibb[0], ibb[1],
					fplot_red, f);
		}
	}

	// render hud
	if (e->head_up_display)
		dump_hud(f);

	f->changed = 1;
}


// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	static double ox = 0, oy = 0;

	if (e->show_lines && m == 0)
	{
		e->dip_xy[0] = x;
		e->dip_xy[1] = y;
		e->dip_pq[2] = 1;
		e->dip_idx = window_to_image(e->dip_pq, e, x, y);
		if (e->dip_idx == 0) // point is on the left, line on the right
			matrix_times_vector(e->dip_abc, e->fmtr, e->dip_pq);
		if (e->dip_idx == 1) // point is on the right, line on the left
			matrix_times_vector(e->dip_abc, e->fm, e->dip_pq);
		f->changed = true;
	} else e->dip_abc[0] = NAN;

	if (m == FTR_BUTTON_LEFT)
		action_offset_viewport(f, x - ox, y - oy);

	if (m == FTR_MASK_SHIFT)
		action_center_contrast_at_point(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	int xy[] = {x, y};
	int idx = subwindow(e, xy, xy+1);

	fprintf(stderr, "button b=%d m=%d\n", b, m);
	if (b == FTR_BUTTON_UP && (m==FTR_MASK_SHIFT || m==FTR_MASK_CONTROL)) {
		action_contrast_span(f, idx, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && ((m==FTR_MASK_SHIFT)||m==FTR_MASK_CONTROL)){
		action_contrast_span(f, idx, 1.3); return; }
	if (b == FTR_BUTTON_DOWN) action_decrease_octave(f, x, y);
	if (b == FTR_BUTTON_UP  ) action_increase_octave(f, x, y);
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
	if (b == FTR_BUTTON_LEFT) e->scroll_domain = idx;
}

void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isprint(k)?k:' ', m, x, y);

	if (k == '+' || k == '=') action_decrease_octave(f, 3*f->w/4, f->h/2);
	if (k == '-') action_increase_octave(f, 3*f->w/4, f->h/2);

	if (k == 'u') action_toggle_hud(f);
	if (k == 'l') action_toggle_lock(f);

	if (k == ',') action_save_shot(f);

	if (k == 'n') action_qauto(f);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// arrows move the viewport
	if (k > 1000) {
		int d[2] = {0, 0};
		int inc = -10;
		if (m & FTR_MASK_SHIFT  ) inc /= 10;
		if (m & FTR_MASK_CONTROL) inc *= 10;
		switch (k) {
		case FTR_KEY_LEFT : d[0] -= inc; break;
		case FTR_KEY_RIGHT: d[0] += inc; break;
		case FTR_KEY_UP   : d[1] -= inc; break;
		case FTR_KEY_DOWN : d[1] += inc; break;
		}
		if (k == FTR_KEY_PAGE_UP)   d[1] = +f->h/3;
		if (k == FTR_KEY_PAGE_DOWN) d[1] = -f->h/3;
		action_offset_viewport(f, d[0], d[1]);
	}
}

#define BAD_BOUND(a,x,b) ((a)<(x)?((x)<(b)?(x):(b)):(a))

void pan_resize(struct FTR *f, int k, int m, int x, int y)
{
	struct pan_state *e = f->userdata;
	e->win_half = f->w / 2;
}

#include "parsenumbers.c"
int main_pan(int c, char *v[])
{

	// process input arguments
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s matrix left.png right.png\n", *v);
		//                          0 1      2        3
		return 1;
	}
	char *matrix_str = v[1];
	char *filename_x = v[2];
	char *filename_y = v[3];

	// read images and fundamental matrix
	struct pan_state e[1];
	e->x[0] = iio_read_image_float_vec(filename_x, 0+e->w, 0+e->h, 0+e->pd);
	e->x[1] = iio_read_image_float_vec(filename_y, 1+e->w, 1+e->h, 1+e->pd);
	read_n_doubles_from_string(e->fm, matrix_str, 9);
	matrix_transpose(e->fmtr, e->fm);
	for (int i = 0; i < 9; i++)
		fprintf(stderr, "fm[%d] = %g\n", i, e->fm[i]);

	// init state
	e->scroll_domain = 0;
	e->dip_abc[0] = NAN;
	e->head_up_display = false;
	e->lock_transform = false;
	e->font = *xfont_9x15;
	e->font = reformat_font(e->font, UNPACKED);

	// open window
	int win_width  = BAD_BOUND(200, e->w[0] + e->w[1], 1024);
	int win_height = BAD_BOUND(200, e->h[0], 800);
	e->win_half = win_width / 2;
	struct FTR f = ftr_new_window(win_width, win_height);

	// bind state to window, initialize handlers, run loop
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	ftr_set_handler(&f, "key"   , pan_key_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "resize", pan_resize);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	return r;
}

int main(int c, char *v[])
{
	return main_pan(c, v);
}
