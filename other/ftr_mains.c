#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "ftr.h"

#include "iio.h"

// read an image into 32-bit BGRA (the "native" X format)
static unsigned char *read_image_uint8_argb(char *fname, int *w, int *h)
{
	int pd;
	unsigned char *x = iio_read_image_uint8_vec(fname, w, h, &pd);
	if (pd == 4) return x;
	unsigned char *y = malloc(4**w**h);
	for (int i = 0; i < *w**h; i++) {
		y[4*i+3] = 255;
		switch(pd) {
		case 1:
			y[4*i+0] = y[4*i+1] = y[4*i+2] = x[i];
			break;
		case 2:
			y[4*i+0] = x[2*i+1];
			y[4*i+1] = y[4*i+2] = y[4*i+3] = x[2*i+0];
			break;
		case 3:
			y[4*i+0] = x[3*i+2];
			y[4*i+1] = x[3*i+1];
			y[4*i+2] = x[3*i+0];
			y[4*i+3] = 255;
			//for (int l = 0; l < 3; l++)
			//	y[4*i+(3-l)] = x[pd*i+l];
			break;
		default:
			for (int l = 0; l < 3; l++)
				y[4*i+l] = x[pd*i+l];
			break;
		}
	}
	free(x);
	return y;
}

static void write_imnage_uint8_argb(char *fname, unsigned char *x, int w, int h)
{
	unsigned char *t = malloc(3*w*h);
	for (int i = 0; i < w*h; i++) {
		t[3*i+0] = x[4*i+2];
		t[3*i+1] = x[4*i+1];
		t[3*i+2] = x[4*i+0];
	}
	iio_save_image_uint8_vec(fname, t, w, h, 3);
	free(t);
}


int main_viewimage(int c, char *v[])
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	int w, h;
	unsigned char *x = read_image_uint8_argb(filename_in, &w, &h);

	struct FTR f = ftr_new_window_with_image_uint8_argb(x, w, h);
	//ftr_change_title(&f, filename_in);
	int r = ftr_loop_run(&f);
	free(x);
	return r;
}

//int main_viewimage2(int c, char *v[])
//{
//	if (c != 2 && c != 1) {
//		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
//		//                          0  1
//		return 1;
//	}
//	char *filename_in = c > 1 ? v[1] : "-";
//
//	int w, h;
//	unsigned char *x = read_image_uint8_argb(filename_in, &w, &h);
//
//	struct FTR f = ftr_new_window(x, w, h);
//
//	int r = ftr_loop_run(&f);
//	free(x);
//	return r;
//}

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b

static void do_inline_crop(unsigned char *x, int *w, int *h, int crop[4])
{
	int from[2], to[2];
	from[0] = BAD_MIN(crop[0], crop[2]);
	from[1] = BAD_MIN(crop[1], crop[3]);
	to[0] = BAD_MAX(crop[0], crop[2]);
	to[1] = BAD_MAX(crop[1], crop[3]);
	int ow = to[0] - from[0]; assert(ow > 0); assert(ow <= *w);
	int oh = to[1] - from[1]; assert(oh > 0); assert(oh <= *h);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		int o_idx = j*ow + i;
		int i_idx = (j+from[1])* *w + i+from[0];
		for (int l = 0; l < 4; l++)
			x[4*o_idx+l] = x[4*i_idx+l];
	}
	*w = ow;
	*h = oh;
}

int main_icrop(int c, char *v[])
{
	if (c != 2 && c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s [in [out]]\n", *v);
		//                          0  1   2
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";
	char *filename_out = c > 2 ? v[2] : "-";

	int w, h;
	unsigned char *x = read_image_uint8_argb(filename_in, &w, &h);

	struct FTR f = ftr_new_window_with_image_uint8_argb(x, w, h);

	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param+0, crop_param+1, NULL, NULL);
	ftr_wait_for_mouse_click(&f, crop_param+2, crop_param+3, NULL, NULL);

	do_inline_crop(x, &w, &h, crop_param);

	write_imnage_uint8_argb(filename_out, x, w, h);
	free(x);
	return 0;
}

int main_pclick(int c, char *v[])
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [in] > pos.txt\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	int w, h;
	unsigned char *x = read_image_uint8_argb(filename_in, &w, &h);

	struct FTR f = ftr_new_window_with_image_uint8_argb(x, w, h);

	int pos[2];
	ftr_wait_for_mouse_click(&f, pos+0, pos+1, NULL, NULL);
	ftr_close(&f);

	printf("%d %d\n", pos[0], pos[1]);

	free(x);
	return 0;
}

static void draw_random(struct FTR *f, int x, int y, int k, int m)
{
	for (int i = 0; i < f->w * f->h * 4; i++)
		f->argb[i] = rand()%256;
	f->changed = 1;
}

int main_fire(int c, char *v[])
{
	int w = 800;
	int h = 600;
	unsigned char *x = malloc(4*w*h);

	struct FTR f = ftr_new_window_with_image_uint8_argb(x, w, h);
	ftr_set_handler(&f, "idle", draw_random);
	ftr_loop_run(&f);

	ftr_close(&f);
	free(x);
	return 0;
}

int main(int c, char *v[])
{
	int (*f)(int,char*[]);
	if (c < 2) return fprintf(stderr, "name a main\n");
	else if (0 == strcmp(v[1], "viewimage")) f = main_viewimage;
	else if (0 == strcmp(v[1], "icrop"))     f = main_icrop;
	else if (0 == strcmp(v[1], "pclick"))    f = main_pclick;
	else if (0 == strcmp(v[1], "fire"))      f = main_fire;
	//else if (0 == strcmp(*v, "simplest"))  f = main_simplest;
	//else if (0 == strcmp(*v, "simplest2")) f = main_simplest2;
	else return fprintf(stderr, "bad main \"%s\"\n", v[1]);
	return f(c-1, v+1);
}
