#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "iio.h"
#include "ftr.h"

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b


static void do_inline_crop_rgb(uint8_t *x, int *w, int *h, int c[4])
{
	int from[2] = {BAD_MIN(c[0], c[2]), BAD_MIN(c[1], c[3])};
	int to[2]   = {BAD_MAX(c[0], c[2]), BAD_MAX(c[1], c[3])};
	int ow = to[0] - from[0]; assert(ow > 0); assert(ow <= *w);
	int oh = to[1] - from[1]; assert(oh > 0); assert(oh <= *h);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		int o_idx = j*ow + i;
		int i_idx = (j+from[1])* *w + i+from[0];
		for (int l = 0; l < 3; l++)
			x[3*o_idx+l] = x[3*i_idx+l];
	}
	*w = ow;
	*h = oh;
}

struct icrop2_state {
	// input data
	uint8_t *original_image;
	int w, h;

	// visualization data
	int octave, offset_x, offset_y;

	// interaction data
	int step;   // step of the interaction (0/1 = before/after first click)
	int px, py; // location of the pointer when in step 0
};

static int inbetween(int a, int b, int x)
{
	return (a <= x && x <= b) || (b <= x && x <= a);
}

void icrop2_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct icrop2_state *e = f->userdata;
	uint8_t *o = e->original_image;

	if (e->step == 0) {
		e->px = x;
		e->py = y;
	}

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < 3; l++)
	{
		int idx = 3*(j * f->w + i) + l;
		int inside = (e->step == 0) ?  (i>=x && j >= y) :
			inbetween(e->px, x, i) && inbetween(e->py, y, j);
		f->rgb[idx] = inside ? o[idx] : o[idx]/2;
	}

	f->changed = 1;
}

// interactive crop, just fancier
int main_icrop2_simple(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	// read image
	int w, h;
	uint8_t (*xx)[3] = iio_read_image_uint8_rgb(filename_in, &w, &h);
	uint8_t *x = *xx;

	// show image in window
	int win[2] = { w, h};//320, 200 };
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, win[0], win[1]);
	struct icrop2_state e = {x,w,h, 0,0,0, 0,0,0};
	f.userdata = &e;

	// set handlers
	ftr_set_handler(&f, "motion", icrop2_motion);

	// read the two corners of the crop rectangle
	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param + 0, crop_param + 1);
	e.step = 1;
	ftr_wait_for_mouse_click(&f, crop_param + 2, crop_param + 3);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	iio_write_image_uint8_vec(filename_out, x, w, h, 3);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

int main_icrop2_panable(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	// read image
	int w, h;
	uint8_t (*xx)[3] = iio_read_image_uint8_rgb(filename_in, &w, &h);
	uint8_t *x = *xx;

	// show image in window
	int win[2] = { w, h};//320, 200 };
	//int win[2] = { 800, 600 };
	struct FTR f = ftr_new_window_with_image_uint8_rgb(x, win[0], win[1]);
	struct icrop2_state e = {x,w,h, 0,0,0, 0,0,0};
	f.userdata = &e;

	// set handlers
	ftr_set_handler(&f, "motion", icrop2_motion);

	// read the two corners of the crop rectangle
	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param + 0, crop_param + 1);
	e.step = 1;
	ftr_wait_for_mouse_click(&f, crop_param + 2, crop_param + 3);

	// perform the crop on the image data
	do_inline_crop_rgb(x, &w, &h, crop_param);

	// write outpuf file
	iio_write_image_uint8_vec(filename_out, x, w, h, 3);

	// cleanup and exit (optional)
	ftr_close(&f);
	free(x);
	return 0;
}

int main_icrop(int c, char *v[]) { return main_icrop2_panable(c, v); }
int main(int c, char *v[]) { return main_icrop(c, v); }
