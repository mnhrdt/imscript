#include "fancy_image.h"

static void fancy_crop(char *fname_out, char *fname_in,
		int x0, int y0, int w, int h)
{
	// open input image
	struct fancy_image *a = fancy_image_open(fname_in, "r");

	// read information from input image
	int tw = 0, th = 0, fmt = 0, bps = 0;
	fancy_image_leak_tiff_info(&tw, &th, &fmt, &bps, a);

	// adapt tile size
	if (tw == a->w) tw = w;
	if (th == a->h) th = h;

	// create output image of the appropriate size and options
	struct fancy_image *b = fancy_image_create(fname_out,
			"w=%d,h=%d,pd=%d,bps=%d,fmt=%d,tw=%d,th=%d",
			w, h, a->pd, bps, fmt, tw, th);

	// fill-in the cropped image
	for (int j = 0; j < b->h; j++)
	for (int i = 0; i < b->w; i++)
	for (int l = 0; l < b->pd; l++)
	{
		float v = fancy_image_getsample(a, x0 + i, y0 + j, l);
		fancy_image_setsample(b, i, j, l, v);
	}

	// close both images
	fancy_image_close(b);
	fancy_image_close(a);
}

#include <stdio.h>
#include <stdlib.h>
int main_fancy_crop(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s x0 y0 w h in.tiff out.tiff\n", *v);
		//        0 1  2  3 4 5       6
		return 1;
	}
	int x0 = atoi(v[1]);
	int y0 = atoi(v[2]);
	int w = atoi(v[3]);
	int h = atoi(v[4]);
	char *filename_in = v[5];
	char *filename_out = v[6];

	fancy_crop(filename_out, filename_in, x0, y0, w, h);

	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_fancy_crop(c, v); }
#endif
