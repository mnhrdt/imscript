// all the flipping functions below have the same signature

static void flip_leftright(
		float *X,    // output image data
		int *WH,     // output image (width,height) in (WH[0],WH[1])
		float *x,    // input image data
		int w,       // input image width
		int h        // input image height
		)
{
	int W = WH[0] = w;
	int H = WH[1] = h;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[j*w + w - 1 - i];
}

static void flip_topdown(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = w;
	int H = WH[1] = h;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[(h-1-j)*w + i];
}

static void flip_transpose(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = h;
	int H = WH[1] = w;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[i*w+j];
}

static void flip_r90(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = h;
	int H = WH[1] = w;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[i*w+(w-1-j)];
}

static void flip_r270(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = h;
	int H = WH[1] = w;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[(h-1-i)*w+j];
}

static void flip_r180(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = w;
	int H = WH[1] = h;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[(h-1-j)*w+(w-1-i)];
}

static void flip_posetrans(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = h;
	int H = WH[1] = w;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[(h-1-i)*w+(w-1-j)];
}

static void flip_identity(float *X, int *WH, float *x, int w, int h)
{
	int W = WH[0] = w;
	int H = WH[1] = h;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
		X[j*W+i] = x[j*w + i];
}

// note: all 8 elements of the D8 are covered now


static char *help_string_name     = "imflip";
static char *help_string_version  = "imflip 1.0\n\nWritten by eml";
static char *help_string_oneliner = "flip an image";
static char *help_string_usage    = "usage:\n\t"
"flip {leftright|topdown|transpose} [in [out]]";
static char *help_string_long     =
"Imflip flips an image in the specified direction.\n"
"\n"
"Usage: imflip OPERATION in.tiff out.png\n"
"   or: imflip OPERATION in.tiff > out.pnm\n"
"   or: cat in.tiff | imflip OPERATION > out.pnm\n"
"\n"
"Operations:\n"
" leftright\t\t\n"
" topdown\t\t\n"
" transpose\t\t\n"
" posetrans\t\t\n"
" r90\t\t\n"
" r180\t\t\n"
" r270\t\t\n"
" identity\t\t\n"
"\n"
"Examples:\n"
" imflip leftright in.png out.png     Flip an image horizontally.\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "help_stuff.c" // functions that print the strings named above
#include "iio.h"
int main_imflip(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// get positional arguments
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t"
			"%s {leftright|topdown|transpose} [in [out]]\n", *v);
		//        0 1                              2   3
		return 1;
	}
	char *op = v[1]; if (*op) *op = tolower(*op);
	char *filename_in = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// read input image and alloc space for output image
	int w, h, pd, wh[2] = {0, 0};
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = malloc(w*h*pd*sizeof*y);

	// flip each channel as requested
	for (int l = 0; l < pd; l++)
	{
		float *X = x + l*w*h;
		float *Y = y + l*w*h;
		if (!strcmp(op, "1"))         flip_identity(Y, wh, X, w, h);
		if (!strcmp(op, "identity"))  flip_identity(Y, wh, X, w, h);
		if (!strcmp(op, "x"))         flip_leftright(Y, wh, X, w, h);
		if (!strcmp(op, "y"))         flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "r"))         flip_r90      (Y, wh, X, w, h);
		if (!strcmp(op, "rr"))        flip_r180     (Y, wh, X, w, h);
		if (!strcmp(op, "rrr"))       flip_r270     (Y, wh, X, w, h);
		if (!strcmp(op, "t"))         flip_transpose(Y, wh, X, w, h);
		if (!strcmp(op, "z"))         flip_posetrans(Y, wh, X, w, h);
		if (!strcmp(op, "leftright")) flip_leftright(Y, wh, X, w, h);
		if (!strcmp(op, "rightleft")) flip_leftright(Y, wh, X, w, h);
		if (!strcmp(op, "lr"))        flip_leftright(Y, wh, X, w, h);
		if (!strcmp(op, "rl"))        flip_leftright(Y, wh, X, w, h);
		if (!strcmp(op, "topdown"))   flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "topbottom")) flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "bottomup"))  flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "updown"))    flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "ud"))        flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "td"))        flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "tb"))        flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "bu"))        flip_topdown  (Y, wh, X, w, h);
		if (!strcmp(op, "transpose")) flip_transpose(Y, wh, X, w, h);
		if (!strcmp(op, "trans"))     flip_transpose(Y, wh, X, w, h);
		if (!strcmp(op, "posetrans")) flip_posetrans(Y, wh, X, w, h);
		if (!strcmp(op, "r90"))       flip_r90      (Y, wh, X, w, h);
		if (!strcmp(op, "r-90"))      flip_r270     (Y, wh, X, w, h);
		if (!strcmp(op, "rm90"))      flip_r270     (Y, wh, X, w, h);
		if (!strcmp(op, "r270"))      flip_r270     (Y, wh, X, w, h);
		if (!strcmp(op, "r180"))      flip_r180     (Y, wh, X, w, h);
	}

	// save and exit
	iio_write_image_float_split(filename_out, y, wh[0], wh[1], pd);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_imflip(c, v); }
#endif
