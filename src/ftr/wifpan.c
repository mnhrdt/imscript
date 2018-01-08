// cc -Ofast wifpan.c iio.o -o wifpan -lX11 -ltiff -ljpeg -lpng -lz -lm
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "iio.h"

#ifndef FTR_BACKEND
#define FTR_BACKEND 'x'
#endif
#include "ftr.h"

#include "xmalloc.c"

#define WHEEL_FACTOR 2
#define MAX_PYRAMID_LEVELS 30

// image file input/output (wrapper around iio) {{{1
static float *read_image_float_rgb(char *fname, int *w, int *h)
{
	int pd;
	float *x = iio_read_image_float_vec(fname, w, h, &pd);
	if (pd == 3) return x;
	float *y = xmalloc(3**w**h*sizeof*y);
	for (int i = 0; i < *w**h; i++) {
		switch(pd) {
		case 1:
			y[3*i+0] = y[3*i+1] = y[3*i+2] = x[i];
			break;
		case 2:
			y[3*i+0] = x[2*i+0];            // R=x
			y[3*i+1] = y[3*i+2] = x[2*i+1]; // G=B=y
			break;
		default:
			assert(pd > 3);
			for (int l = 0; l < 3; l++)
				y[3*i+l] = x[pd*i+l];
			break;
		}
	}
	free(x);
	return y;
}

// data structure for the image viewer
// this data goes into the "userdata" field of the FTR window structure
struct pan_state {
	// 1. image data
	int w, h;
	float *frgb;

	// 2. view port parameters
	double zoom_factor, offset_x, offset_y;
	double a, b;
	double aaa[3], bbb[3];

	// 3. image pyramid
	float *pyr_rgb[MAX_PYRAMID_LEVELS];
	int pyr_w[MAX_PYRAMID_LEVELS], pyr_h[MAX_PYRAMID_LEVELS];

	// 4. roi
	int roi; // 0=nothing 1=dftwindow 2=rawfourier 3=ppsmooth
	int roi_x, roi_y, roi_w;
};

// change of coordinates: from window "int" pixels to image "double" point
static void window_to_image(double p[2], struct pan_state *e, int i, int j)
{
	p[0] = e->offset_x + i / e->zoom_factor;
	p[1] = e->offset_y + j / e->zoom_factor;
}

// change of coordinates: from image "double" point to window "int" pixel
static void image_to_window(int i[2], struct pan_state *e, double x, double y)
{
	i[0] = floor(x * e->zoom_factor - e->offset_x);
	i[1] = floor(y * e->zoom_factor - e->offset_y);
}

static float getsample_0(float *x, int w, int h, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[3*(j*w+i)+l];
}

// evaluate a frgb image at a point
static void interpolate_at(float *out, float *x, int w, int h, float p, float q)
{
	int i = p;
	int j = q;
	if (i < 0 || i >= w) { out[0] = out[1] = out[2] = 0; return; }
	if (j < 0 || j >= h) { out[0] = out[1] = out[2] = 0; return; }
	out[0] = x[3*(j*w+i)+0];
	out[1] = x[3*(j*w+i)+1];
	out[2] = x[3*(j*w+i)+2];
	//return x[3*(j*w+i)+l];
	//out[0] = getsample_0(x, w, h, (int)p, (int)q, 0);
	//out[1] = getsample_0(x, w, h, (int)p, (int)q, 1);
	//out[2] = getsample_0(x, w, h, (int)p, (int)q, 2);
}

// evaluate the shadow of a frgb image (the g and b channels are ignored)
static
void interpolate_at_s(float *out, float *x, int w, int h, float p, float q)
{
	float oo = getsample_0(x, w, h, 0+(int)p, 0+(int)q, 0);
	float ox = getsample_0(x, w, h, 1+(int)p, 0+(int)q, 0);
	float oy = getsample_0(x, w, h, 0+(int)p, 1+(int)q, 0);
	float od = getsample_0(x, w, h, 1+(int)p, 1+(int)q, 0);
	float out_s = NAN;
	float out_si = 0;
	if (isfinite(oo)) {
		out_s = 0;
		if (isfinite(ox) && isfinite(oy) && isfinite(od)) {
			float gx = ox - oo + od - oy;
			float gy = oy - oo + od - ox;
			//float nx[3] = {1, 0, gx};
			//float ny[3] = {0, 1, gy};
			float nn[3] = {-gx/2, -gy/2, 1};
			float ss[3] = {-1, -1, -1};
			out_s = nn[0]*ss[0] + nn[1]*ss[1] + nn[2]*ss[2];
		}
		out_si = oo > 300 ? 255 : 0;
		out_s = 127 + 10 * out_s;
	}
	out[0] = out[1] = out[2] = out_s;
	out[0] = (out[1]+out_si)/2;
	out[2] = (out[2]+255-out_si)/2;
}

// evaluate the value a position (p,q) in image coordinates
static void pixel(float *out, struct pan_state *e, double p, double q)
{
	if (e->zoom_factor > 0.9999)
		interpolate_at(out, e->frgb, e->w, e->h, p, q);
	else {
		//static int first_run = 1;
		//if (first_run) {
		//	fprintf(stderr, "create pyramid\n");
		//	void create_pyramid(struct pan_state *e);
		//	create_pyramid(e);
		//	first_run = 0;
		//}
		if(p<0||q<0){out[0]=out[1]=out[2]=0;return;}
		int s = round(log2(1/e->zoom_factor)) - 1;
		if (s < 0) s = 0;
		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
		int sfac = 1<<(s+1);
		int w = e->pyr_w[s];
		int h = e->pyr_h[s];
		float *rgb = e->pyr_rgb[s];
		interpolate_at(out, rgb, w, h, p/sfac, q/sfac);
	}
}

// pixel shadow
//static void pixel_s(float *out, struct pan_state *e, double p, double q)
//{
//	if (e->zoom_factor > 0.9999)
//		interpolate_at_s(out, e->frgb, e->w, e->h, p, q);
//	else {
//		//static int first_run = 1;
//		//if (first_run) {
//		//	fprintf(stderr, "create pyramid\n");
//		//	void create_pyramid(struct pan_state *e);
//		//	create_pyramid(e);
//		//	first_run = 0;
//		//}
//		if(p<0||q<0){out[0]=out[1]=out[2]=0;return;}
//		int s = -0 - log(e->zoom_factor) / log(2);
//		if (s < 0) s = 0;
//		if (s >= MAX_PYRAMID_LEVELS) s = MAX_PYRAMID_LEVELS-1;
//		int sfac = 1<<(s+1);
//		int w = e->pyr_w[s];
//		int h = e->pyr_h[s];
//		float *rgb = e->pyr_rgb[s];
//		interpolate_at_s(out, rgb, w, h, p/sfac, q/sfac);
//	}
//}

static void action_print_value_under_cursor(struct FTR *f, int x, int y)
{
	if (x<f->w && x>=0 && y<f->h && y>=0) {
		struct pan_state *e = f->userdata;
		double p[2];
		window_to_image(p, e, x, y);
		float c[3];
		interpolate_at(c, e->frgb, e->w, e->h, p[0], p[1]);
		printf("%g\t%g\t: %g\t%g\t%g\n", p[0], p[1], c[0], c[1], c[2]);
	}
}

static void action_offset_viewport(struct FTR *f, int dx, int dy)
{
	struct pan_state *e = f->userdata;
	e->offset_x -= dx/e->zoom_factor;
	e->offset_y -= dy/e->zoom_factor;

	f->changed = 1;
}

static void action_reset_zoom_and_position(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	e->zoom_factor = 1;
	e->offset_x = 0;
	e->offset_y = 0;
	e->a = 1;
	e->b = 0;
	e->bbb[0] = e->bbb[1] = e->bbb[2] = 0;

	e->roi = 0;
	e->roi_x = f->w / 2;
	e->roi_y = f->h / 2;
	e->roi_w = 73; // must be odd

	f->changed = 1;
}

static void action_contrast_change(struct FTR *f, float afac, float bshift)
{
	struct pan_state *e = f->userdata;

	e->a *= afac;
	e->b += bshift;

	f->changed = 1;
}

static void action_qauto(struct FTR *f)
{
	struct pan_state *e = f->userdata;

	float m = INFINITY, M = -m;
	int pid = 3;
	for (int i = 0; i < 3 * e->pyr_w[pid] * e->pyr_h[pid]; i++)
	{
		float g = e->pyr_rgb[pid][i];
		m = fmin(m, g);
		M = fmax(M, g);
	}

	e->a = 255 / ( M - m );
	e->b = 255 * m / ( m - M );
	e->bbb[0] = e->b;
	e->bbb[1] = e->b;
	e->bbb[2] = e->b;

	f->changed = 1;
}

static void action_toggle_roi(struct FTR *f, int x, int y, int dir)
{
	struct pan_state *e = f->userdata;
	e->roi = (e->roi + (dir?-1:1)) % 4;
	fprintf(stderr, "ROI SWITCH(%d) = %d\n", dir, e->roi);
	e->roi_x = x - e->roi_w / 2;
	e->roi_y = y - e->roi_w / 2;
	f->changed = 1;
}

static void action_roi_embiggen(struct FTR *f, int d)
{
	struct pan_state *e = f->userdata;
	e->roi_w += d;
	fprintf(stderr, "ROI %d\n", e->roi_w);
	f->changed = 1;
}


static void action_move_roi(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;
	if (e->roi)
	{
		e->roi_x = x - e->roi_w / 2;
		e->roi_y = y - e->roi_w / 2;
		f->changed = 1;
	}
}

static void action_center_contrast_at_point(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	double p[2];
	window_to_image(p, e, x, y);
	float c[3];
	pixel(c, e, p[0], p[1]);
	float C = (c[0] + c[1] + c[2])/3;

	e->bbb[0] = 127.5 - e->a * c[0];
	e->bbb[1] = 127.5 - e->a * c[1];
	e->bbb[2] = 127.5 - e->a * c[2];

	e->b = 127.5 - e->a * C;

	f->changed = 1;
}

static void action_contrast_span(struct FTR *f, float factor)
{
	struct pan_state *e = f->userdata;

	float c = (127.5 - e->b)/ e->a;
	float ccc[3];
	for(int l=0;l<3;l++) ccc[l] = (127.5 - e->bbb[l]) / e->a;
	e->a *= factor;
	e->b = 127.5 - e->a * c;
	for(int l=0;l<3;l++) e->bbb[l] = 127.5 - e->a * ccc[l];

	f->changed = 1;
}

static void action_change_zoom_by_factor(struct FTR *f, int x, int y, double F)
{
	struct pan_state *e = f->userdata;

	double c[2];
	window_to_image(c, e, x, y);

	e->zoom_factor *= F;
	e->offset_x = c[0] - x/e->zoom_factor;
	e->offset_y = c[1] - y/e->zoom_factor;
	fprintf(stderr, "\t zoom changed %g\n", e->zoom_factor);

	f->changed = 1;
}

static void action_reset_zoom_only(struct FTR *f, int x, int y)
{
	struct pan_state *e = f->userdata;

	action_change_zoom_by_factor(f, x, y, 1/e->zoom_factor);
}



static void action_increase_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, WHEEL_FACTOR);
}

static void action_decrease_zoom(struct FTR *f, int x, int y)
{
	action_change_zoom_by_factor(f, x, y, 1.0/WHEEL_FACTOR);
}

static bool insideP(int w, int h, int i, int j)
{
	return i>=0 && j>=0 && i<w && j<h;
}

#define OMIT_BLUR_MAIN
#include "blur.c"

#define OMIT_PPSMOOTH_MAIN
#include "ppsmooth.c"

static void ppsmooth_vec(float *y, float *x, int w, int h, int pd)
{
	float *x_split = malloc(w*h*pd*sizeof*x);
	float *y_split = malloc(w*h*pd*sizeof*x);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		x_split[l*w*h+(j*w+i)] = x[pd*(j*w+i)+l];

	ppsmooth_split(y_split, x_split, w, h, pd);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		y[pd*(j*w+i)+l] = y_split[l*w*h+(j*w+i)];

	free(x_split);
	free(y_split);
}

// note: assumes 3-dimensional pixels
static void transform_roi_buffers_old(float *y, float *x, int n)
{
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
		float *xx = x + 3*(j*n + i);
		float *yy = y + 3*(j*n + i);
		float g = hypot(xx[0], hypot(xx[1],xx[2]))/sqrt(3);
		yy[0] = sqrt(g/255)*255;
		yy[1] = g/2;
		yy[2] = g;
	}
	if (false) return;
	float x_p[3*n*n];
	float x_s[3*n*n];
	ppsmooth_vec(x_p, x, n, n, 3);
	for (int i = 0; i < 3*n*n; i++)
		x_s[i] = x[i] - x_p[i];
	float param[1] = {4};
	blur_2d(y, x_p, n, n, 3, "laplace", param, 1);
	for (int i = 0; i < 3*n*n; i++)
		y[i] += x_s[i];
}

static void transform_roi_buffers(float *y, float *x, int n, int roi)
{
	float x_p0[3*n*n], *x_p = x_p0;
	ppsmooth_vec(x_p0, x, n, n, 3);
	if (roi == 2) x_p = x;
	if (roi == 3) { for (int i = 0; i < 3*n*n; i++) y[i]=x_p0[i]; return; }
	float *c = xmalloc(n*n*sizeof*c);
	float *ys = xmalloc(n*n*sizeof*c);
	fftwf_complex *fc = fftwf_xmalloc(n*n*sizeof*fc);
	for (int l = 0; l < 3; l++)
	{
		for (int i = 0; i < n*n; i++)
			c[i] = x_p[3*i+l];
		fft_2dfloat(fc, c, n, n);
		for (int i = 0; i < n*n; i++)
			ys[i] = 255*(log(cabs(fc[i])/255)+0.5)/5;
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			int ii = (i + n/2) % n;
			int jj = (j + n/2) % n;
			y[3*(j*n+i)+l] = ys[jj*n+ii];
		}
	}
	free(ys);
	free(c);
	fftwf_free(fc);
	//float param[1] = {0.5};
	//blur_2d(y, y, n, n, 3, "cauchy", param, 1);
}


// dump the image acording to the state of the viewport
static void pan_exposer(struct FTR *f, int b, int m, int x, int y)
{
	struct pan_state *e = f->userdata;

	// expose the whole image
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2];
		window_to_image(p, e, i, j);
		float c[3];
		pixel(c, e, p[0], p[1]);
		//if (!i && !j)
		//	fprintf(stderr, "first pixel (%g %g){%g %g %g}\n",
		//			p[0], p[1], c[0], c[1], c[2]);
		unsigned char *cc = f->rgb + 3 * (j * f->w + i);
		for (int l = 0; l < 3; l++)
		{
			if (!isfinite(c[l]))
				cc[l] = 0;
			else {
				//float g = e->a * c[l] + e->b;
				float g = e->a * c[l] + e->bbb[l];
				if      (g < 0)   cc[l] = 0  ;
				else if (g > 255) cc[l] = 255;
				else              cc[l] = g  ;
			}
		}
	}

	// if ROI, expose the roi
	if (e->roi)
	{
		float buf_in [3 * e->roi_w * e->roi_w];
		float buf_out[3 * e->roi_w * e->roi_w];
		for (int j = 0; j < e->roi_w; j++)
		for (int i = 0; i < e->roi_w; i++)
		{
			// i,h       = position inside ROI
			// ii,jj     = position inside viewport
			// p[0],p[1] = position inside image domain
			int ii = i + e->roi_x - e->roi_w/2;
			int jj = j + e->roi_y - e->roi_w/2;
			double p[2];
			window_to_image(p, e, ii, jj);
			float *c = buf_in + 3 * (j * e->roi_w + i);
			pixel(c, e, p[0], p[1]);
		}
		transform_roi_buffers(buf_out, buf_in, e->roi_w, e->roi);
		for (int j = 0; j < e->roi_w; j++)
		for (int i = 0; i < e->roi_w; i++)
		{
			int ii = i + e->roi_x - e->roi_w/2;
			int jj = j + e->roi_y - e->roi_w/2;
			float *c = buf_out + 3 * (j * e->roi_w + i);
			if (insideP(f->w, f->h, ii, jj))
			{
				uint8_t *cc = f->rgb + 3 * (jj * f->w + ii);
				for (int l = 0; l < 3; l++)
				{
					float g = e->a * c[l] + e->bbb[l];
					if      (g < 0)   cc[l] = 0  ;
					else if (g > 255) cc[l] = 255;
					else              cc[l] = g  ;
					//cc[l] = c[l];
				}
			}
		}
	}

	// mark shit as changed
	f->changed = 1;
}

// update offset variables by dragging
static void pan_motion_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "motion b=%d m=%d (%d %d)\n", b, m, x, y);

	static double ox = 0, oy = 0;

	if (m & FTR_BUTTON_LEFT)   action_offset_viewport(f, x - ox, y - oy);
	if (m & FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (m & FTR_MASK_SHIFT)    action_center_contrast_at_point(f, x, y);

	action_move_roi(f, x, y);

	ox = x;
	oy = y;
}

static void pan_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	//fprintf(stderr, "button b=%d m=%d\n", b, m);
	struct pan_state *e = f->userdata;

	if (e->roi && b == FTR_BUTTON_UP) { action_roi_embiggen(f,+2); return; }
	if (e->roi && b == FTR_BUTTON_DOWN){action_roi_embiggen(f,-2); return; }
	if (b == FTR_BUTTON_UP && m & FTR_MASK_SHIFT) {
		action_contrast_span(f, 1/1.3); return; }
	if (b == FTR_BUTTON_DOWN && m & FTR_MASK_SHIFT) {
		action_contrast_span(f, 1.3); return; }
	if (b == FTR_BUTTON_RIGHT && m & FTR_MASK_CONTROL) {
		action_reset_zoom_only(f, x, y); return; }
	if (b == FTR_BUTTON_MIDDLE) action_print_value_under_cursor(f, x, y);
	if (b == FTR_BUTTON_DOWN)   action_increase_zoom(f, x, y);
	if (b == FTR_BUTTON_UP  )   action_decrease_zoom(f, x, y);
	if (b == FTR_BUTTON_RIGHT)  action_reset_zoom_and_position(f);
}

static void key_handler_print(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "key pressed %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
}

static void pan_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	//fprintf(stderr, "PAN_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
	//		k, isalpha(k)?k:' ', m, x, y);

	//if (k == '+') action_increase_zoom(f, f->w/2, f->h/2);
	//if (k == '-') action_decrease_zoom(f, f->w/2, f->h/2);
	if (k == '+') action_change_zoom_by_factor(f, f->w/2, f->h/2, 2);
	if (k == '-') action_change_zoom_by_factor(f, f->w/2, f->h/2, 0.5);
	if (k == 'p') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.1);
	if (k == 'm') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.1);
	if (k == 'P') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1.006);
	if (k == 'M') action_change_zoom_by_factor(f, f->w/2, f->h/2, 1/1.006);

	//if (k == 'a') action_contrast_change(f, 1.3, 0);
	//if (k == 'A') action_contrast_change(f, 1/1.3, 0);
	//if (k == 'b') action_contrast_change(f, 1, 1);
	//if (k == 'B') action_contrast_change(f, 1, -1);
	if (k == 'n') action_qauto(f);
	if (k == 'r') action_toggle_roi(f, x, y, m&FTR_MASK_SHIFT);

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

	// if 'k', do weird things
	if (k == 'k') {
		fprintf(stderr, "setting key_handler_print\n");
		ftr_set_handler(f, "key", key_handler_print);
	}
}

static void zoom_out_by_factor_two_avg(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < 3; l++)
	{
		float a[4], ax = 0;
		a[0] = getsample_0(in, iw, ih, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, 2*i+1, 2*j+1, l);
		int cx = 0;
		for (int i = 0; i < 4; i++)
			if (!isnan(a[i])) {
				ax += a[i];
				cx += 1;
			}
		out[3*(ow*j + i)+l] = cx ? ax/cx : NAN;//(a[0] + a[1] + a[2] + a[3])/4;
	}
}

static void zoom_out_by_factor_two_max(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < 3; l++)
	{
		float a[4];
		a[0] = getsample_0(in, iw, ih, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, 2*i+1, 2*j+1, l);
		out[3*(ow*j + i)+l] = fmax(fmax(a[0],a[1]),fmax(a[2],a[3]));
	}
}

// type of a "zoom-out" function
typedef void (*zoom_out_function_t)(float*,int,int,float*,int,int);

#include "smapa.h"
SMART_PARAMETER_SILENT(ZOMAX,0)
SMART_PARAMETER_SILENT(PYRSAVE,0)

static void create_pyramid(struct pan_state *e)
{
	zoom_out_function_t z = zoom_out_by_factor_two_avg;
	if (ZOMAX() > 0)
		z = zoom_out_by_factor_two_max;
	for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
	{
		int      lw   = s ? e->pyr_w  [s-1] : e->w   ;
		int      lh   = s ? e->pyr_h  [s-1] : e->h   ;
		float   *lrgb = s ? e->pyr_rgb[s-1] : e->frgb;
		int      sw   = ceil(lw / 2.0);
		int      sh   = ceil(lh / 2.0);
		float   *srgb = xmalloc(3 * sw * sh * sizeof*srgb);
		z(srgb, sw, sh, lrgb, lw, lh);
		e->pyr_w[s]   = sw;
		e->pyr_h[s]   = sh;
		e->pyr_rgb[s] = srgb;
	}
	if (PYRSAVE() > 0) {
		for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
		{
			char buf[FILENAME_MAX];
			snprintf(buf, FILENAME_MAX, "/tmp/pyrsave_%d.tiff", s);
			iio_write_image_float_vec(buf, e->pyr_rgb[s],
					e->pyr_w[s], e->pyr_h[s], 3);
			fprintf(stderr, "saved image \"%s\"\n", buf);
		}
	}
}

static void free_pyramid(struct pan_state *e)
{
	for (int s = 0; s < MAX_PYRAMID_LEVELS; s++)
		free(e->pyr_rgb[s]);
}


#define BAD_MIN(a,b) a<=b?a:b
int main_wifpan(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	// read image
	struct pan_state e[1];
	e->frgb = read_image_float_rgb(filename_in, &e->w, &e->h);
	create_pyramid(e);

	// open window
	struct FTR f = ftr_new_window(BAD_MIN(e->w,1000), BAD_MIN(e->h,800));
	f.userdata = e;
	action_reset_zoom_and_position(&f);
	ftr_set_handler(&f, "expose", pan_exposer);
	ftr_set_handler(&f, "motion", pan_motion_handler);
	ftr_set_handler(&f, "button", pan_button_handler);
	ftr_set_handler(&f, "key"   , pan_key_handler);
	int r = ftr_loop_run(&f);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(e->frgb);
	free_pyramid(e);
	return r - 1;
}

int main(int c, char *v[])
{
	return main_wifpan(c, v);
}
