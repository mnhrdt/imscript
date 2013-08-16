#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xmalloc.c"
#include "xfopen.c"
#include "fail.c"


struct bitmap_font {
	int number_of_glyphs;
	int width;
	int height;
	enum {UNPACKED, PACKED, ZRLE, HUFFMAN, LZW, RLEHUFFMAN} packing;
	unsigned char *data;
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

// PACKED => UNPACKED
static struct bitmap_font unpack_font(struct bitmap_font *fp)
{
	struct bitmap_font fu = *fp;
	int usize = fu.width * fu.height * fu.number_of_glyphs;
	int psize = usize / 8;
	fu.data = unpack_bit_data(fp->data, psize);
	assert(fp->packing == PACKED);
	fu.packing = UNPACKED;
	return fu;
}

// UNPACKED => PACKED
static struct bitmap_font pack_font(struct bitmap_font *fu)
{
	struct bitmap_font fp = *fu;
	int usize = fp.width * fp.height * fp.number_of_glyphs;
	int psize = usize / 8;
	fp.data = pack_bit_data(fu->data, usize);
	assert(fu->packing == UNPACKED);
	fp.packing = PACKED;
	return fp;
}

static void dump_font_as_parseable_c_struct(char *filename,
		struct bitmap_font *ufont, char *name)
{
	struct bitmap_font font = pack_font(ufont);
	FILE *f = xfopen(filename, "w");
	int datasize = font.number_of_glyphs * font.width * font.height;
	fprintf(f, "unsigned char %s_data[] = {\n\t", name);
	for (int i = 0; i < datasize/8; i++)
		fprintf(f, "%d,%s", (unsigned char)font.data[i],
				(i&&!(i%16))?"\n\t":" ");
	fprintf(f, "0\n};\n");
	fprintf(f, "struct bitmap_font %s[1] = "
			"{{%d, %d, %d, PACKED, %s_data}};\n",
		name, font.number_of_glyphs, font.width, font.height, name);
	xfclose(f);
	free(font.data);
}

static int get_font_bit(struct bitmap_font *f, int c, int i, int j)
{
	assert(f->packing == UNPACKED);
	if (c < 0 || c > f->number_of_glyphs)
		return 0;
	return f->data[(c*f->height + j)*f->width + i];
}

static void set_font_bit(struct bitmap_font *f, int c, int i, int j)
{
	assert(f->packing == UNPACKED);
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

// scan a line of a BDF file
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

// artificial and arbitrary limit to the number of glyphs
#define NGLYPHS 256

static void font_fill_from_bdf(struct bitmap_font *font, char *fname)
{
	//fprintf(stderr, "going to parse BDF file \"%s\"\n", fname);
	FILE *f = xfopen(fname, "r");
	{
		float bdfversion;
		int r = fscanf(f, "STARTFONT %g\n", &bdfversion);
		//fprintf(stderr, "BDF version \"%g\"\n", bdfversion);
	}
	font->data = NULL;
	font->number_of_glyphs = NGLYPHS;
	font->width = font->height = -1;
	int bmp_line = -2, glyph = -1;
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
			size_t ds = num[0] * num[1] * NGLYPHS;
			font->data = xmalloc(ds);
			font->packing = UNPACKED;
			memset(font->data, 0, ds);
		}
		//if (0 == strcmp(tag, "STARTCHAR")) fprintf(stderr, "BDF \"%s\"\n", line);
		if (0 == strcmp(tag, "ENCODING") && nn == 1) glyph = num[0];
		if (0 == strcmp(tag, "BITMAP")) bmp_line = -1;
		if (0 == strcmp(tag, "ENDCHAR")) bmp_line = -2;
		if (bmp_line >= 0)
		{
			//fprintf(stderr, "glyph %d got bitmap line %d:\t\"%s\"", glyph, bitmap_line, tag);
			char bits[1000];
			hex_string_to_bits(bits, tag);
			//for (int i = 0; i < font->width; i++)
			//	fprintf(stderr, bits[i]?"M":"_");
			//fprintf(stderr, " \n");
			for (int i = 0; i < font->width; i++)
				if (bits[i])
					set_font_bit(font, glyph, i, bmp_line);
		}
		if (bmp_line > -2) bmp_line += 1;
	}
	xfclose(f);
	if (!font->data)
		fail("failed to parse BDF file \"%s\"\n", fname);
}

static void put_pixel(float *x, int w, int h, int pd, int i, int j, float *c)
{
	if (j*w + i < w*h)
		for (int l = 0; l < pd; l++)
			x[(w*j+i)*pd+l] = c[l];
}

static void put_string_in_float_image(float *x, int w, int h, int pd,
		int posx, int posy, float *color, int kerning,
		struct bitmap_font *font, char *string)
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

// fontu cdump name {packed|unpacked|zrle|huffman|lzw|mhuffman} [in.bdf [out.c]]
// fontu puts [-f bdf] [-c color] "string" [in.png [out.png]]

static int main_cdump(int c, char **v)
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s name {packing} [in.bdf [out.c]]\n", *v);
		//        0  1    2         3       4
		return 1;
	}
	char *name = v[1];
	char *packing = v[2];
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	struct bitmap_font f;
	font_fill_from_bdf(&f, filename_in);
	dump_font_as_parseable_c_struct(filename_out, &f, name);

	return 0;
}

// fontu puts [-f bdf] [-c color] "string" [in.png [out.png]]

#include "iio.h"

static int main_puts(int c, char **v)
{
	if (c != 5 && c != 6 && c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s font.bdf px py \"string\" [in [out]]\n", *v);
		//        0  1       2  3    4         5   6
		return 1;
	}
	char *bdf = v[1];
	int posx = atoi(v[2]);
	int posy = atoi(v[3]);
	char *text = v[4];
	char *filename_in = c > 5 ? v[5] : "-";
	char *filename_out = c > 6 ? v[6] : "-";

	struct bitmap_font f;
	font_fill_from_bdf(&f, bdf);

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	float black[10] = {0};
	put_string_in_float_image(x, w, h, pd, posx, posy, black, 0, &f, text);

	iio_save_image_float_vec(filename_out, x, w, h, pd);

	free(f.data);
	free(x);
	return 0;
}

int main(int c, char **v)
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "cdump")) return main_cdump(c-1, v+1);
	else if (0 == strcmp(v[1], "puts")) return main_puts(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s {cdump|puts} params\n", *v);
	       return 1;
	}
}
