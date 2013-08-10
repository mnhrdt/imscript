#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xmalloc.c"
#include "xfopen.c"
#include "fail.c"


// Why do we need two structures for fonts?
// 1. packed is space-efficient to embed into a C program
// 2. unpacked is very easy to use inside the program
//
// Note: they are actually the same structure but the interperation is different

struct unpacked_bitmap_font {
	int number_of_glyphs;
	int width;
	int height;
	unsigned char *data; // each char is either 0 or 1
};

struct packed_bitmap_font {
	int number_of_glyphs;
	int width;
	int height;
	unsigned char *data; // each char codifies 8 bits, without padding
};


#define SETBIT(x,i) ((x)|=(1<<(i)))
#define GETBIT(x,i) ((x)&(1<<(i)))

// u: array of chars containing boolean values
static unsigned char *pack_bit_data(unsigned char *u, int nu)
{
	int np = nu / 8;
	if (np*8 != nu)
		fail("can not unpack an odd number of bits");
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

// p: array of chars to be separated into bits
static unsigned char *unpack_bit_data(unsigned char *p, int np)
{
	int nu = 8*np;
	unsigned char *u = xmalloc(nu);
	for (int i = 0; i < np; i++)
		for (int j = 0; j < 8; j++)
			u[8*i+j] = GETBIT(p[i], j);
	return u;
}

static struct unpacked_bitmap_font unpack_font(struct packed_bitmap_font *fp)
{
	struct unpacked_bitmap_font fu; memcpy(&fu, fp, sizeof*fp);
	int usize = fu.width * fu.height * fu.number_of_glyphs;
	int psize = usize / 8;
	fu.data = unpack_bit_data(fp->data, psize);
	return fu;
}

static struct packed_bitmap_font pack_font(struct unpacked_bitmap_font *fu)
{
	struct packed_bitmap_font fp; memcpy(&fp, fu, sizeof*fu);
	int usize = fp.width * fp.height * fp.number_of_glyphs;
	int psize = usize / 8;
	fp.data = pack_bit_data(fu->data, usize);
	return fp;
}

static void dump_font_as_parseable_c_struct(char *filename,
		struct unpacked_bitmap_font *ufont, char *name)
{
	struct packed_bitmap_font font = pack_font(ufont);
	FILE *f = xfopen(filename, "w");
	fprintf(f, "struct packed_bitmap_font %s = { %d, %d, %d, {\n", name,
			font.number_of_glyphs, font.width, font.height);
	int datasize = font.number_of_glyphs * font.width * font.height;
	for (int i = 0; i < datasize/8; i++)
		fprintf(f, "%d,\n", (unsigned char)font.data[i]);
	fprintf(f, "0}\n};\n");
	xfclose(f);
	free(font.data);
}

static int get_font_bit(struct unpacked_bitmap_font *f, int c, int i, int j)
{
	if (c < 0 || c > f->number_of_glyphs)
		return 0;
	return f->data[(c*f->height + j)*f->width + i];
}

static void set_font_bit(struct unpacked_bitmap_font *f, int c, int i, int j)
{
	assert(f->data);
	assert(i >= 0);
	assert(j >= 0);
	assert(i < f->width);
	assert(j < f->height);
	if (c >= 0 && c < f->number_of_glyphs
			&& i >= 0 && j >= 0
			&& i < f->width && j < f->height)
		f->data[(c*f->height + j)*f->width + i] = 1;
}

static int get_tagged_numbers(char *tag, int *nums, char *line)
{
	while (*line && !isspace(*line))
		*tag++ = *line++;
	return sscanf(line, " %d %d %d %d\n", nums, nums+1, nums+2, nums+2);
}

#define EVENP(x) (!((x)&1))

// always fills-in a whole set of bytes, whose quantity is returned
static int hex_string_to_bits(char *bits, char *string)
{
	int nt = strlen(string);
	assert(EVENP(nt));
	int cx = 0;
	for (int j = 0; j < nt; j += 2)
	{
		int c1 = string[j];
		int c2 = string[j+1];
		int nc1 = isdigit(c1) ? c1 - '0' : c1 - 'A' + 10;
		int nc2 = isdigit(c2) ? c2 - '0' : c2 - 'A' + 10;
		int nc = 16*nc1 + nc2;
		//fprintf(stderr, "<%02x>", nc);
		for (int i = 0; i < 8; i++)
		{
			bits[cx++] = GETBIT(nc, 7-i);
			//if (GETBIT(nc,7-i))
			//	fprintf(stderr, "X");
			//else
			//	fprintf(stderr, ".");
		}
	}
	return nt/2;
}

// max line length on a bdf file
#define BDFLL 4096

static void font_fill_from_bdf(struct unpacked_bitmap_font *font, char *fname)
{
	//fprintf(stderr, "goint to parse BDF file \"%s\"\n", fname);
	FILE *f = xfopen(fname, "r");
	{
		float bdfversion;
		int r = fscanf(f, "STARTFONT %g\n", &bdfversion);
		//fprintf(stderr, "BDF version \"%g\"\n", bdfversion);
	}
	font->data = NULL;
	font->number_of_glyphs = 256;
	font->width = font->height = -1;
	int bitmap_line = -2, glyph = -1;
	while (1) {
		char line[BDFLL], tag[BDFLL]={0}, *sl = fgets(line, BDFLL, f);
		if (!sl) break;
		int num[4], nn = get_tagged_numbers(tag, num, line);
		////if (nn >= 0)
		//{
		//	fprintf(stderr, "BDF tagged_numbers(%d) = %s", nn, tag);
		//	for (int i = 0; i < nn; i++)
		//		fprintf(stderr, " %d", num[i]);
		//	fprintf(stderr, "\n");
		//}
		if (0 == strcmp(tag, "FONTBOUNDINGBOX") && nn == 4)
		{
			if (font->width > 0)
				fail("BDF bounding box defined twice");
			font->width = num[0];
			font->height = num[1];
			size_t ds = font->width * font->height * font->number_of_glyphs;
			font->data = xmalloc(ds);
			memset(font->data, 0, ds);
		}
		//if (0 == strcmp(tag, "STARTCHAR")) fprintf(stderr, "BDF \"%s\"\n", line);
		if (0 == strcmp(tag, "ENCODING") && nn == 1) glyph = num[0];
		if (0 == strcmp(tag, "BITMAP")) bitmap_line = -1;
		if (0 == strcmp(tag, "ENDCHAR")) bitmap_line = -2;
		if (bitmap_line >= 0)
		{
			//fprintf(stderr, "glyph %d got bitmap line %d:\t\"%s\"", glyph, bitmap_line, tag);
			char bits[1000];
			hex_string_to_bits(bits, tag);
			//for (int i = 0; i < font->width; i++)
			//	fprintf(stderr, bits[i]?"M":"_");
			//fprintf(stderr, " \n");
			for (int i = 0; i < font->width; i++)
				if (bits[i])
					set_font_bit(font, glyph, i, bitmap_line);
		}
		if (bitmap_line > -2) bitmap_line += 1;
	}
	xfclose(f);
}

static void put_pixel(float *x, int w, int h, int pd, int i, int j, float *c)
{
	if (j*w + i < w*h)
		for (int l = 0; l < pd; l++)
			x[(w*j+i)*pd+l] = c[l];
}

static void put_string_in_float_image(float *x, int w, int h, int pd,
		int posx, int posy, float *color, int kerning,
		struct unpacked_bitmap_font *font, char *string)
{
	while (1)
	{
		int c = *string++;
		if (!c) break;
		if (c > 0 && c < font->number_of_glyphs)
		{
			//fprintf(stderr, "putting glyph \"%d\" '%c'\n", c, c);
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

#include "font_6x12.c"

static float *put_watermark(float *x, int w, int h, int pd, char *url,
		int *ow, int *oh)
{
	int vmargin = 94;
	*ow = w;
	*oh = h + vmargin;
	size_t sf = sizeof(float);
	float *y = xmalloc(*ow * *oh * pd * sf);
	for (int i = 0; i < *ow * *oh * pd; i++)
		y[i] = 255;
	//memset(y, 0xff, *ow * *oh * pd * sf);
	memcpy(y, x, w*h*pd*sf);

	struct packed_bitmap_font pf[1];
	pf->number_of_glyphs = 256;
	pf->width = 6;
	pf->height = 12;
	pf->data = font_data_6x12;
	struct unpacked_bitmap_font f = unpack_font(pf);
	//font_fill_from_bdf(f, "/home/coco/.fonts2/6x12.bdf");
	//dump_font_as_parseable_c_struct("/tmp/6x12.c", f, "6x12");
	float color[pd];
	for (int i = 0; i < pd; i++) color[i] = 0;
	//color[pd-1] = 255;
	put_string_in_float_image(y, *ow, *oh, pd, 10, h+20, color, 0, &f, url);
	return y;
}

static float *remove_watermark(float *x, int w, int h, int pd, int *ow, int *oh)
{
	fail("watermark removal not implemented");
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
		y = put_watermark(x, w, h, pd, url, &ow, &oh);
	else if (operation[0] == 'r')
		y = remove_watermark(x, w, h, pd, &ow, &oh);
	else goto exiterr;

	iio_save_image_float_vec(filename_out, y, ow, oh, pd);

	free(x);
	free(y);
	return 0;
}
