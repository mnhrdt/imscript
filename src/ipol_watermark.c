#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// utility function that returns a valid pointer to memory
static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new) {
		fprintf(stderr, "ERROR: out of memory when requesting "
			       "%zu bytes\n", size);
		abort();
	}
	return new;
}

// data structure for storing a bitmap font
struct bitmap_font {
	int number_of_glyphs;
	int width;
	int height;
	enum {UNPACKED, PACKED, ZRLE} packing;
	unsigned char *data;
};

// macros for accessing the individual bits of an integer
#define SETBIT(x,i) ((x)|=(1<<(i)))
#define GETBIT(x,i) ((x)&(1<<(i)))

// pack an array of boolean (char) values into indivitual bits
static unsigned char *pack_bit_data(unsigned char *u, int nu)
{
	int np = nu / 8;
	if (np*8 != nu) {
		fprintf(stderr, "can not unpack an odd number of bools\n");
		abort();
	}
	unsigned char *p = xmalloc(np);
	for (int i = 0; i < np; i++)
	{
		p[i] = 0;
		for (int j = 0; j < 8; j++)
			if (u[8*i+j])
				SETBIT(p[i], j);
	}
	return p;
}

// unpack a bit field into an array of boolean (char) values
static unsigned char *unpack_bit_data(unsigned char *p, int np)
{
	int nu = 8*np;
	unsigned char *u = xmalloc(nu);
	for (int i = 0; i < np; i++)
		for (int j = 0; j < 8; j++)
			u[8*i+j] = GETBIT(p[i], j);
	return u;
}

// unpack font data
static struct bitmap_font unpack_font(struct bitmap_font *fp)
{
	struct bitmap_font fu; memcpy(&fu, fp, sizeof*fp);
	int usize = fu.width * fu.height * fu.number_of_glyphs;
	int psize = usize / 8;
	fu.data = unpack_bit_data(fp->data, psize);
	assert(fp->packing == PACKED);
	fu.packing = UNPACKED;
	return fu;
}

// pack font data
static struct bitmap_font pack_font(struct bitmap_font *fu)
{
	struct bitmap_font fp; memcpy(&fp, fu, sizeof*fu);
	int usize = fp.width * fp.height * fp.number_of_glyphs;
	int psize = usize / 8;
	fp.data = pack_bit_data(fu->data, usize);
	assert(fu->packing == UNPACKED);
	fp.packing = PACKED;
	return fp;
}


// fixed font 6x12 from X windows
static unsigned char font_data_6x12[] = {
	0, 0, 84, 64, 4, 68, 64, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 4, 65, 16, 0, 1, 0, 0, 160,
	40, 10, 0, 0, 0, 0, 0, 0, 0, 0, 202, 167, 40, 159, 2, 0, 0, 64, 56, 85,
	225, 80, 149, 67, 0, 0, 0, 76, 19, 66, 8, 89, 6, 0, 0, 0, 8, 69, 33,
	84, 137, 5, 0, 0, 64, 16, 4, 0, 0, 0, 0, 0, 0, 128, 16, 132, 32, 8, 4,
	129, 0, 0, 32, 16, 4, 130, 32, 4, 33, 0, 0, 0, 16, 149, 67, 56, 21, 1,
	0, 0, 0, 0, 4, 241, 17, 4, 0, 0, 0, 0, 0, 0, 0, 0, 134, 49, 0, 0, 0, 0,
	0, 240, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 134, 1, 0, 0, 0, 64, 8, 66, 8,
	66, 0, 0, 0, 0, 48, 146, 36, 73, 18, 3, 0, 0, 0, 16, 6, 65, 16, 132, 3,
	0, 0, 0, 56, 17, 132, 16, 194, 7, 0, 0, 0, 124, 16, 194, 64, 145, 3, 0,
	0, 0, 32, 140, 146, 124, 8, 2, 0, 0, 0, 124, 193, 3, 65, 145, 3, 0, 0,
	0, 48, 66, 240, 68, 145, 3, 0, 0, 0, 124, 16, 130, 16, 4, 1, 0, 0, 0,
	56, 81, 228, 68, 145, 3, 0, 0, 0, 56, 81, 228, 65, 136, 1, 0, 0, 0, 0,
	128, 97, 0, 134, 1, 0, 0, 0, 0, 128, 97, 0, 134, 49, 0, 0, 0, 0, 8, 33,
	16, 8, 0, 0, 0, 0, 0, 192, 7, 124, 0, 0, 0, 0, 0, 0, 2, 129, 16, 2, 0,
	0, 0, 0, 56, 17, 66, 16, 0, 1, 0, 0, 0, 56, 81, 87, 117, 129, 3, 0, 0,
	0, 56, 81, 244, 69, 81, 4, 0, 0, 0, 60, 146, 228, 72, 210, 3, 0, 0, 0,
	56, 81, 16, 4, 145, 3, 0, 0, 0, 60, 146, 36, 73, 210, 3, 0, 0, 0, 124,
	65, 240, 4, 193, 7, 0, 0, 0, 124, 65, 240, 4, 65, 0, 0, 0, 0, 56, 81,
	16, 100, 145, 3, 0, 0, 0, 68, 81, 244, 69, 81, 4, 0, 0, 0, 56, 4, 65,
	16, 132, 3, 0, 0, 0, 112, 8, 130, 32, 137, 1, 0, 0, 0, 68, 73, 49, 20,
	73, 4, 0, 0, 0, 4, 65, 16, 4, 193, 7, 0, 0, 0, 68, 91, 21, 69, 81, 4,
	0, 0, 0, 68, 209, 84, 101, 81, 4, 0, 0, 0, 56, 81, 20, 69, 145, 3, 0,
	0, 0, 60, 81, 244, 4, 65, 0, 0, 0, 0, 56, 81, 20, 85, 137, 5, 0, 0, 0,
	60, 81, 244, 20, 73, 4, 0, 0, 0, 56, 81, 224, 64, 145, 3, 0, 0, 0, 124,
	4, 65, 16, 4, 1, 0, 0, 0, 68, 81, 20, 69, 145, 3, 0, 0, 0, 68, 81, 20,
	41, 10, 1, 0, 0, 0, 68, 81, 20, 85, 149, 2, 0, 0, 0, 68, 145, 66, 40,
	81, 4, 0, 0, 0, 68, 145, 66, 16, 4, 1, 0, 0, 0, 124, 16, 66, 8, 193, 7,
	0, 0, 224, 8, 130, 32, 8, 130, 224, 0, 0, 0, 4, 130, 64, 32, 8, 4, 0,
	0, 224, 32, 8, 130, 32, 8, 226, 0, 0, 64, 40, 17, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 124, 0, 32, 16, 8, 0, 0, 0, 0, 0, 0, 0, 0, 128, 3,
	121, 145, 7, 0, 0, 0, 4, 193, 19, 69, 209, 3, 0, 0, 0, 0, 128, 19, 5,
	145, 3, 0, 0, 0, 64, 144, 23, 69, 145, 7, 0, 0, 0, 0, 128, 19, 61, 129,
	3, 0, 0, 0, 48, 146, 112, 8, 130, 0, 0, 0, 0, 0, 128, 19, 69, 145, 7,
	57, 0, 0, 4, 193, 19, 69, 81, 4, 0, 0, 0, 16, 128, 65, 16, 132, 3, 0,
	0, 0, 64, 0, 6, 65, 16, 36, 49, 0, 0, 4, 65, 148, 28, 73, 4, 0, 0, 0,
	24, 4, 65, 16, 132, 3, 0, 0, 0, 0, 192, 82, 85, 85, 5, 0, 0, 0, 0, 64,
	51, 69, 81, 4, 0, 0, 0, 0, 128, 19, 69, 145, 3, 0, 0, 0, 0, 192, 19,
	69, 209, 19, 4, 0, 0, 0, 128, 23, 69, 145, 7, 65, 0, 0, 0, 64, 51, 5,
	65, 0, 0, 0, 0, 0, 128, 23, 56, 208, 3, 0, 0, 0, 16, 196, 71, 16, 4, 6,
	0, 0, 0, 0, 64, 20, 69, 153, 5, 0, 0, 0, 0, 64, 20, 69, 10, 1, 0, 0, 0,
	0, 64, 20, 85, 149, 2, 0, 0, 0, 0, 64, 164, 16, 74, 4, 0, 0, 0, 0, 64,
	20, 69, 10, 33, 4, 0, 0, 0, 192, 135, 16, 194, 7, 0, 0, 128, 16, 4, 33,
	16, 4, 129, 0, 0, 64, 16, 4, 65, 16, 4, 65, 0, 0, 32, 16, 4, 129, 16,
	4, 33, 0, 0, 0, 0, 128, 84, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 65, 16, 4, 1, 0, 0, 0, 0, 132, 83,
	21, 149, 67, 0, 0, 0, 48, 146, 112, 8, 82, 3, 0, 0, 0, 0, 64, 165, 68,
	74, 5, 0, 0, 0, 68, 202, 71, 124, 4, 1, 0, 0, 0, 16, 4, 1, 16, 4, 1, 0,
	0, 192, 9, 140, 36, 49, 144, 3, 0, 0, 160, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	120, 97, 90, 166, 161, 7, 0, 0, 163, 48, 128, 3, 0, 0, 0, 0, 0, 0, 0,
	0, 165, 20, 10, 5, 0, 0, 0, 0, 0, 240, 65, 16, 0, 0, 0, 0, 0, 0, 224,
	0, 0, 0, 0, 0, 0, 120, 97, 91, 150, 161, 7, 0, 0, 240, 1, 0, 0, 0, 0,
	0, 0, 0, 35, 73, 12, 0, 0, 0, 0, 0, 0, 0, 16, 196, 71, 16, 192, 7, 0,
	132, 130, 16, 14, 0, 0, 0, 0, 0, 6, 66, 32, 6, 0, 0, 0, 0, 0, 0, 128,
	16, 2, 0, 0, 0, 0, 0, 0, 0, 0, 64, 20, 69, 217, 21, 4, 0, 224, 93, 215,
	101, 81, 20, 5, 0, 0, 0, 0, 0, 195, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	128, 24, 132, 65, 16, 14, 0, 0, 0, 0, 0, 132, 66, 0, 14, 0, 0, 0, 0, 0,
	0, 0, 0, 64, 161, 80, 74, 1, 0, 194, 32, 8, 10, 163, 120, 8, 2, 0, 194,
	32, 8, 10, 5, 33, 4, 7, 0, 3, 33, 16, 11, 163, 120, 8, 2, 0, 0, 0, 16,
	0, 65, 8, 145, 3, 0, 2, 1, 56, 81, 244, 69, 81, 4, 0, 8, 1, 56, 81,
	244, 69, 81, 4, 0, 132, 2, 56, 81, 244, 69, 81, 4, 0, 86, 3, 56, 81,
	244, 69, 81, 4, 0, 128, 2, 56, 81, 244, 69, 81, 4, 0, 132, 66, 56, 81,
	244, 69, 81, 4, 0, 0, 0, 120, 69, 241, 20, 69, 7, 0, 0, 0, 56, 81, 16,
	4, 145, 131, 24, 2, 1, 124, 65, 240, 4, 193, 7, 0, 8, 1, 124, 65, 240,
	4, 193, 7, 0, 132, 2, 124, 65, 240, 4, 193, 7, 0, 128, 2, 124, 65, 240,
	4, 193, 7, 0, 2, 1, 56, 4, 65, 16, 132, 3, 0, 8, 1, 56, 4, 65, 16, 132,
	3, 0, 132, 2, 56, 4, 65, 16, 132, 3, 0, 128, 2, 56, 4, 65, 16, 132, 3,
	0, 0, 0, 56, 146, 116, 73, 146, 3, 0, 86, 3, 68, 209, 84, 101, 81, 4,
	0, 2, 1, 56, 81, 20, 69, 145, 3, 0, 8, 1, 56, 81, 20, 69, 145, 3, 0,
	132, 2, 56, 81, 20, 69, 145, 3, 0, 86, 3, 56, 81, 20, 69, 145, 3, 0,
	128, 2, 56, 81, 20, 69, 145, 3, 0, 0, 0, 0, 145, 66, 40, 17, 0, 0, 0,
	0, 57, 89, 85, 85, 147, 19, 0, 2, 1, 68, 81, 20, 69, 145, 3, 0, 8, 1,
	68, 81, 20, 69, 145, 3, 0, 132, 2, 68, 81, 20, 69, 145, 3, 0, 128, 2,
	68, 81, 20, 69, 145, 3, 0, 8, 1, 68, 145, 66, 16, 4, 1, 0, 0, 0, 8,
	142, 36, 73, 142, 0, 0, 0, 0, 56, 81, 82, 36, 81, 3, 0, 0, 32, 16, 128,
	3, 121, 145, 7, 0, 0, 128, 16, 128, 3, 121, 145, 7, 0, 0, 64, 40, 128,
	3, 121, 145, 7, 0, 0, 96, 53, 128, 3, 121, 145, 7, 0, 0, 0, 40, 128, 3,
	121, 145, 7, 0, 0, 64, 40, 132, 3, 121, 145, 7, 0, 0, 0, 0, 128, 67,
	57, 133, 7, 0, 0, 0, 0, 128, 19, 5, 145, 131, 24, 0, 32, 16, 128, 19,
	61, 129, 3, 0, 0, 128, 16, 128, 19, 61, 129, 3, 0, 0, 64, 40, 128, 19,
	61, 129, 3, 0, 0, 0, 40, 128, 19, 61, 129, 3, 0, 0, 32, 16, 128, 65,
	16, 132, 3, 0, 0, 128, 16, 128, 65, 16, 132, 3, 0, 0, 64, 40, 128, 65,
	16, 132, 3, 0, 0, 0, 40, 128, 65, 16, 132, 3, 0, 128, 66, 40, 144, 23,
	69, 145, 3, 0, 0, 96, 53, 64, 51, 69, 81, 4, 0, 0, 32, 16, 128, 19, 69,
	145, 3, 0, 0, 128, 16, 128, 19, 69, 145, 3, 0, 0, 64, 40, 128, 19, 69,
	145, 3, 0, 0, 96, 53, 128, 19, 69, 145, 3, 0, 0, 0, 40, 128, 19, 69,
	145, 3, 0, 0, 0, 0, 4, 240, 1, 4, 0, 0, 0, 0, 0, 128, 151, 85, 211, 3,
	0, 0, 32, 16, 64, 20, 69, 145, 3, 0, 0, 128, 16, 64, 20, 69, 145, 3, 0,
	0, 64, 40, 64, 20, 69, 145, 3, 0, 0, 0, 40, 64, 20, 69, 145, 3, 0, 0,
	128, 16, 64, 20, 69, 10, 33, 4, 0, 0, 4, 193, 19, 69, 209, 19, 4, 0, 0,
	40, 64, 20, 69, 10, 33, 4, 0
};
static struct bitmap_font font_6x12[1] = {{256, 6, 12, PACKED, font_data_6x12}};


// get the pixel (i,j) of glyph c from the font f
static int get_font_bit(struct bitmap_font *f, int c, int i, int j)
{
	if (c < 0 || c > f->number_of_glyphs)
		return 0;
	return f->data[(c*f->height + j)*f->width + i];
}


// set the pixel (i,j) of glyph c of the font f
static void set_font_bit(struct bitmap_font *f, int c, int i, int j)
{
	if (c >= 0 && c < f->number_of_glyphs
			&& i >= 0 && i < f->width
			&& j >= 0 && j < f->height)
		f->data[(c*f->height + j)*f->width + i] = 1;
}

// set pixel (i,j) to color c
static void put_pixel(float *x, int w, int h, int pd, int i, int j, float *c)
{
	if (j*w + i < w*h)
		for (int l = 0; l < pd; l++)
			x[(w*j+i)*pd+l] = c[l];
}

// print a string into an image
static void put_string_in_float_image(float *x, int w, int h, int pd,
		int posx, int posy, float *color, int kerning,
		struct bitmap_font *font, char *string)
{
	int n = strlen(string);
	for (int k = 0; k < n; k++)
	{
		int c = string[k];
		if (c >= 0 && c < font->number_of_glyphs)
		{
			for (int i = 0; i < font->width; i++)
			for (int j = 0; j < font->height; j++)
				if (get_font_bit(font, c, i, j))
				{
					int ii = posx + i;
					int jj = posy + j;
					put_pixel(x, w, h, pd, ii, jj, color);
				}
		}
		posx += font->width + kerning;
	}
}


#define MIN_WIDTH 432
#define UMARGIN 84

#define LINEH 16
#define LINE_1 "This image was created by an IPOL demo"
#define LINE_2 url
#define LINE_3 "Please cite the corresponding article if you use this image."
#define LINE_4 "To remove this comment, upload the image to http://dev.ipol.im/unmark"

#define MAGIC_PX 0
#define MAGIC_PY 0
#define MAGIC_UINT16 57089
#define MAGIC_BPP 1

static void unfold_byte(int *buf, int byte)
{
	assert(byte >= 0);
	assert(byte < 256);
	for (int i = 0; i < 8; i++)
		buf[i] = (bool)GETBIT(byte, i);
}

static int fold_byte(int *buf)
{
	int r = 0;
	for (int i = 0; i < 8; i++)
		if (buf[i])
			SETBIT(r, i);
	assert(r >= 0);
	assert(r < 256);
	return r;
}

static void get_buf_from_code(int *buf, uint16_t code[5])
{
	for (int i = 0; i < 5; i++)
	{
		//fprintf(stderr, "code[%d] = %hd\n", i, (unsigned)code[i]);
		unfold_byte(buf + 16*i + 0, code[i]%256);
		unfold_byte(buf + 16*i + 8, code[i]/256);
	}
	//for (int i = 0 ; i < 80; i++)
	//	fprintf(stderr, "buf[%d] = %d\n", i, buf[i]);
}

static void get_code_from_buf(uint16_t code[5], int *buf)
{
	//for (int i = 0 ; i < 80; i++)
	//	fprintf(stderr, "buf[%d] = %d\n", i, buf[i]);
	for (int i = 0; i < 5; i++)
	{
		int lo = fold_byte(buf + 16*i + 0);
		int hi = fold_byte(buf + 16*i + 8);
		code[i] = 256*hi + lo;
		//fprintf(stderr, "code[%d] = %hd\n", i, (unsigned)code[i]);
	}
}

// put the IPOL watermark
static float *put_ipol_watermark(float *x, int w, int h, int pd, char *url,
		int *oow, int *ooh)
{
	// 1. the white margin
	int umargin = UMARGIN;
	int rmargin = 0;
	if (w < MIN_WIDTH)
		rmargin = MIN_WIDTH - w;

	int ow = w + rmargin;
	int oh = h + umargin;
	size_t sf = sizeof(float);
	float *y = xmalloc(ow * oh * pd * sf);
	for (int i = 0; i < ow * oh * pd; i++)
		y[i] = 255; // set background white
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		y[(ow*(j+umargin)+i)*pd+l] = x[(w*j+i)*pd+l];

	// 2. the watermark text
	struct bitmap_font f = unpack_font(font_6x12);
	float color[pd];
	for (int i = 0; i < pd; i++) color[i] = 0;

	put_string_in_float_image(y,ow,oh,pd, 5, 5, color, 0, &f,LINE_1);
	if (pd == 3) {color[0] = 0; color[1] = 0; color[2] = 255;} // blue
	put_string_in_float_image(y,ow,oh,pd, 25,5+LINEH, color, 0, &f,LINE_2);
	if (pd == 3) {color[0] = 0; color[1] = 0; color[2] = 0;}   // black
	put_string_in_float_image(y,ow,oh,pd, 5,5+2*LINEH, color, 0, &f,LINE_3);
	if (pd == 3) {color[0] = 0; color[1] = 170; color[2] = 0;} // green
	put_string_in_float_image(y,ow,oh,pd, 5,5+4*LINEH, color, 0, &f,LINE_4);

	// 3. the "magic" code
	uint16_t code[5] = {MAGIC_UINT16, 0, umargin, w, h};
	int buf[80] = {0};
	get_buf_from_code(buf, code);
	for (int k = 0; k < 80; k++)
	{
		int i = MAGIC_PX + k;
		int j = MAGIC_PY;
		y[(j*ow+i)*pd] -= buf[k];
	}

	*oow = ow;
	*ooh = oh;
	return y;
}

// remove the IPOL watermark
static float *remove_ipol_watermark(float *x, int w, int h, int pd,
		int *ow, int *oh)
{
	// 1. find magic code, and associated data
	int buf[80];
	for (int k = 0; k < 80; k++)
	{
		int i = MAGIC_PX + k;
		int j = MAGIC_PY;
		buf[k] = round(255-x[(j*w+i)*pd]);
	}
	uint16_t code[5];
	get_code_from_buf(code, buf);

	if (*code != MAGIC_UINT16)
	{
		fprintf(stderr, "could not find magic unwatermarking code\n");
		abort();
	}

	fprintf(stderr, "orig = %d %d\n", w, h);
	fprintf(stderr, "crop data = %d %d %d %d\n",
			code[1], code[2], code[3], code[4]);

	if (code[1] >= w) {fprintf(stderr,"inconsistent crop\n");abort();}
	if (code[2] >= h) {fprintf(stderr,"ynconsistent crop\n");abort();}
	if (code[1]+code[3]>w){fprintf(stderr,"inconsistent crop\n");abort();}
	if (code[2]+code[4]>h){fprintf(stderr,"ynconsistent krop\n");abort();}

	// 2. perform the actual crop
	*ow = code[3];
	*oh = code[4];
	float *y = xmalloc(*ow * *oh * pd * sizeof*y);
	for (int i = 0; i < code[3]; i++)
	for (int j = 0; j < code[4]; j++)
	for (int l = 0; l < pd; l++)
	{
		int ii = code[1] + i;
		int jj = code[2] + j;
		y[(*ow*j + i)*pd+l] = x[(w*jj + ii)*pd+l];
	}

	return y;
}





#include "iio.h"

int main(int c, char **v)
{
	if (c != 5) {
exiterr:	fprintf(stderr, "usage:\n\t"
				"%s {put|remove} in.png URL out.png\n", *v);
		//                0  1           2      3   4
		return 1;
	}
	char *operation = v[1];
	char *filename_in = v[2];
	char *url = v[3];
	char *filename_out = v[4];

	int w, h, pd, ow, oh;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *y;

	if (operation[0] == 'p')
		y = put_ipol_watermark(x, w, h, pd, url, &ow, &oh);
	else if (operation[0] == 'r')
		y = remove_ipol_watermark(x, w, h, pd, &ow, &oh);
	else goto exiterr;

	iio_save_image_float_vec(filename_out, y, ow, oh, pd);

	free(x);
	free(y);
	return 0;
}
