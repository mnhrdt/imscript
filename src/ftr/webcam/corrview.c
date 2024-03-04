#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "iio.h" // for saving the screenshot
#include "ftr.h"
#include "cam.h"

#include "seconds.c"

#define OMIT_BLUR_MAIN
#include "blur.c"

#define OMIT_PPSMOOTH_MAIN
#include "ppsmooth.c"

#define OMIT_MORSI_MAIN
#include "morsi.c"

struct kam_state {
	struct camera *c;
	uint8_t *prev;

	int diff_mode;
	int simplest_color_balance;

	int corr_mode;
	int corr_w;
};

static void filterview_autocorr(float *y, float *x, int w, int h)
{
	float *z = xmalloc(w*h*sizeof*z);
	memset(y, w*h*sizeof*y, 0);

	// morphological filtering
	int *E = build_disk(5.1);
	morsi_bothat(z, x, w, h, E);

	// laplacian
	for (int j = 1; j < h-1; j++)
	for (int i = 1; i < w-1; i++)
		y[j*w+i] = 127 - 1*(-4*z[j*w+i]
			+z[(j+1)*w+i+1]+z[(j+1)*w+i-1]
			+z[(j-1)*w+i+1]+z[(j-1)*w+i-1] );



	for (int i = 0; i < w*h; i++)
		y[i] = 4*z[i];

	for (int i = 0; i < w*h; i++)
		y[i] = fmax(0, fmin(y[i], 255));

	free(E);
	free(z);
}

static void filterview_autocorr_old(float *y, float *x, int w, int h)
{
	// laplacian
	memset(y, w*h*sizeof*y, 0);
	for (int j = 1; j < h-1; j++)
	for (int i = 1; i < w-1; i++)
		y[j*w+i] = 127 - 1*(-4*x[j*w+i]
			+x[(j+1)*w+i+1]+x[(j+1)*w+i-1]
			+x[(j-1)*w+i+1]+x[(j-1)*w+i-1] );

	// morphological filtering
	float *z = xmalloc(w*h*sizeof*z);
	int *E = build_disk(5.1);
	morsi_enhance(z, y, w, h, E);


	for (int i = 0; i < w*h; i++)
		y[i] = 1*z[i];

	for (int i = 0; i < w*h; i++)
		y[i] = fmax(0, fmin(y[i], 255));

	free(E);
	free(z);
}

static void autocorrelation_inplace(struct kam_state *e, float *x, int n)
{
	//for (int j = 0; j < n; j++)
	//for (int i = 0; i < n; i++)
	//	x[j*n+i] = 255 - x[j*n+i];

	//double tic = seconds();

	float *c = malloc(n*n*sizeof*c);
	ppsmooth_split(c, x, n, n, 1);
	float *ys = malloc(n*n*sizeof*ys);
	fftwf_complex *fc = fftwf_xmalloc(n*n*sizeof*fc);

	fft_2dfloat(fc, c, n, n);
	for (int i = 0; i < n*n; i++)
		fc[i] = cabs(fc[i]);
	ifft_2dfloat(c, fc, n, n);
	for (int i = 0; i < n*n; i++)
		ys[i] = c[i];
	//float fac = n*3;
	for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
		int ii = (i + n/2) % n;
		int jj = (j + n/2) % n;
		//float norm = hypot(i-n/2-1, j-n/2-1) / fac;
		c[j*n+i] = ys[jj*n+ii] * 1;//norm;
	}

	if (e->diff_mode)
		filterview_autocorr(x, c, n, n);
	else
		filterview_autocorr_old(x, c, n, n);

	free(ys);
	free(c);
	fftwf_free(fc);

	//double tac = seconds();
	//fprintf(stderr, "autocorr n=%d ms=%g\n", n, 1000*(tac - tic));
}

static void kam_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	static double tic = 0;
	double tac = seconds();
	fprintf(stderr, "tictac = %g ms (%g fps)\n",1000*(tac-tic),1/(tac-tic));
	tic = tac;

	(void)unused_x; (void)unused_y;
	struct kam_state *e = f->userdata;

	// workspace: dark blue background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 100 * (2 == i%3);

	struct camera *c = e->c;
	camera_grab_rgb(c);

	// leave a symmetric margin
	int ox = (f->w - c->w) / 2;
	int oy = (f->h - c->h) / 2;
	assert(ox >= 0);
	assert(oy >= 0);

	for (int j = 0; j < c->h; j++)
	for (int i = 0; i < c->w; i++)
	for (int k = 0; k < 3; k++)
	{
		float v = c->rgb [3*(j*c->w+i)+k];
		if (v > 255) v = 255;
		if (v < 0) v = 0;
		f->rgb[3*(f->w*(j+oy)+i+ox)+k] = v;
	}


	if (e->corr_mode == 0) // none
		;
	if (e->corr_mode == 1) // gray
	{
		float *aimg = xmalloc(e->corr_w * e->corr_w * sizeof*aimg);
		int i0 = f->w/2 - e->corr_w/2 - ox;
		int j0 = f->h/2 - e->corr_w/2 - oy;
		for (int j = 0; j < e->corr_w; j++)
		for (int i = 0; i < e->corr_w; i++)
		{
			float r = f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+0];
			float g = f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+1];
			float b = f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+2];
			float v = (r+g+b)/3;
			aimg[j*e->corr_w + i] = v;
		}
		autocorrelation_inplace(e, aimg, e->corr_w);
		for (int j = 0; j < e->corr_w; j++)
		for (int i = 0; i < e->corr_w; i++)
		for (int k = 0; k < 3; k++)
			f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+k] =
					aimg[j*e->corr_w+i];
		free(aimg);
	}
	if (e->corr_mode == 2) // separate rgb
	{
		int n = e->corr_w;
		float *aimg = xmalloc(3*n*n*sizeof*aimg);
		int i0 = f->w/2 - n/2 - ox;
		int j0 = f->h/2 - n/2 - oy;
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		for (int k = 0; k < 3; k++)
		{
			float v = f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+k];
			aimg[k*n*n + j*n + i] = v;
		}
		for (int k = 0; k < 3; k++)
			autocorrelation_inplace(e, aimg + k*n*n, n);
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		for (int k = 0; k < 3; k++)
			f->rgb[3*(f->w*(j+j0+oy)+i+i0+ox)+k] =
					aimg[k*n*n + j*n+i];
		free(aimg);
	}


	//if (!e->diff_mode)
	//{
	//	for (int j = 0; j < c->h; j++)
	//		memcpy(
	//				f->rgb+ 3*(f->w*(j+oy)+ox),
	//				c->rgb+ 3*(c->w*j),
	//				3*c->w
	//		      );
	//} else {
	//	for (int j = 0; j < c->h; j++)
	//	for (int i = 0; i < c->w; i++)
	//	for (int k = 0; k < 3   ; k++)
	//	{
	//		float v1 = c->rgb [3*(j*c->w+i)+k];
	//		float v2 = e->prev[3*(j*c->w+i)+k];
	//		float v = 3*(v2 - v1) + 127;
	//		if (v > 255) v = 255;
	//		if (v < 0) v = 0;
	//		f->rgb[3*(f->w*(j+oy)+i+ox)+k] = v;
	//	}
	//	memcpy(e->prev, c->rgb, 3 * c->w * c->h);
	//}
}

static void action_take_jpeg_screenshot(struct kam_state *e)
{
	static int c = 0;
	char n[FILENAME_MAX];
	snprintf(n, FILENAME_MAX, "webcam_%d.jpg", c);
	iio_write_image_uint8_vec(n, e->c->rgb, e->c->w, e->c->h, 3);
	fprintf(stderr, "wrote sreenshot on file \"%s\"\n", n);
	c += 1;
}

void kam_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);
	fprintf(stderr, "KEY k=%d ('%c') m=%d x=%d y=%d\n", k, k, m, x, y);

	struct kam_state *e = f->userdata;

	if (k == 'c') e->corr_mode = (1 + e->corr_mode) % 3;
	if (k == 'd') e->diff_mode = !e->diff_mode;
	if (k == 'b') e->simplest_color_balance = !e->simplest_color_balance;
	if (k == 'j') action_take_jpeg_screenshot(e);

	if  (k == '\033' || k == 'q') // ESC or q
		ftr_notify_the_desire_to_stop_this_loop(f, 1);
}

int main()
{
	// state of the program
	struct kam_state e[1];

	// camera stuff
	int w = 640;//800;
	int h = 480;//600;
	e->c = camera_begin("/dev/video0", w, h); // XXX: depends on cam!

	e->prev = malloc(2*3 * w * h);
	e->diff_mode = 0;

	e->corr_mode = 1;
	e->corr_w = 256;

	// window stuff
	struct FTR f = ftr_new_window(w + 100, h + 100);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , kam_key_handler);
	ftr_set_handler(&f, "idle", kam_exposer);
	int r = ftr_loop_run(&f);
	ftr_close(&f);

	free(e->prev);
	camera_end(e->c);


	return 0;
}
