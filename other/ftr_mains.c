#include "ftr.h"

#include "iio.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

int main_viewimage(int c, char *v[])
{
	if (c != 2 && c != 1) {
		fprintf(stderr, "usage:\n\t%s [image]\n", *v);
		//                          0  1
		return 1;
	}
	char *filename_in = c > 1 ? v[1] : "-";

	int w, h;
	unsigned char *x = iio_read_image_uint8_rgba(filename_in, &w, &h);

	struct FTR f = ftr_new_window_with_image_uint8_rgba(x, w, h);
	ftr_change_title(&f, filename_in);
	int r = ftr_loop_run(&f);
	free(x);
	return r;
}

#define BAD_MIN(a,b) a<b?a:b
#define BAD_MAX(a,b) a>b?a:b

static void do_inline_crop(unsigned char *x, int *w, int *h, int crop[4])
{
	int from[2], to[2];
	from[0] = BAD_MIN(crop[0], crop[2]);
	from[1] = BAD_MIN(crop[1], crop[3]);
	to[0] = BAD_MAX(crop[0], crop[2]);
	to[1] = BAD_MAX(crop[1], crop[3]);
	int ow = to[0] - from[0]; assert(ow > 0); assert(ow <= w);
	int oh = to[1] - from[1]; assert(oh > 0); assert(oh <= h);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		int o_idx = j*ow + i;
		int i_idx = (j+from[1])*w + i+from[0];
		for (int l = 0; l < 4; l++)
			x[o_idx] = x[i_idx];
	}
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
	unsigned char *x = iio_read_image_uint8_rgba(filename_in, &w, &h);

	struct FTR f = ftr_new_window_with_image_uint8_rgba(x, w, h);

	int crop_param[4];
	ftr_wait_for_mouse_click(&f, crop_param+0, crop_param+1, NULL, NULL);
	ftr_wait_for_mouse_click(&f, crop_param+2, crop_param+3, NULL, NULL);

	do_inline_crop(x, &w, &h, crop_param);

	iio_save_image_uint8_vec(filename_out, x, w, h);
	free(x);
	return 0;
}

int main(int c, char *v[])
{
	int (*f)(int,char*[]);
	if (false) ;
	else if (0 == strcmp(*v, "viewimage")) f = main_viewimage;
	else if (0 == strcmp(*v, "icrop"))     f = main_icrop;
	//else if (0 == strcmp(*v, "simplest"))  f = main_simplest;
	//else if (0 == strcmp(*v, "simplest2")) f = main_simplest2;
	else return fprintf(stderr, "bad main\n");
	return f(c-1, v+1);
}
