#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xmalloc.c"
#include "xfopen.c"
#include "fail.c"

enum font_data_format {
	UNPACKED, // array of chars with boolean values
	PACKED,   // array of chars with arbitrary values (with the same bitmap)
	RLE,      // binary run-length encoding
	PCX,      // binary run-length encoding
	RLET,      // binary run-length encoding, transposed
	//HUFFMAN,  // huffman code of "PACKED"
	//RLEH,     // huffman code of "RLE"
	X85,      // x85 encoding of "PACKED"
	X85RLE,   // x85 encoding of "RLE"
	X85RLET,   // x85 encoding of "RLE", transposed
	//X85HUF,   // x85 encoding of "HUFFMAN"
	//X85RLEH,  // x85 encoding of "RLEHUFFMAN" (efficient for C dumps)
	DIFF,
	XOR,
	RLEDIFF,
	RLEXOR,
	RLEPCX,
	RLEXORPCX,
	PCXX85,
	RLEPCXX85,
	RLEXORPCXX85,
};

// a font is a three-dimensional binary image
// its bits can be stored in different formats
struct bitmap_font {
	int number_of_glyphs;
	int width;
	int height;

	enum font_data_format packing;
	int ndata; // when UNPACKED: ndata = number_of_glyphs * width * height
	char *name; // optional field
	unsigned char *data;
};

//#define SETBIT(x,i) ((x)|=(1<<(i)))
//#define GETBIT(x,i) (bool)((x)&(1<<(i)))
//
//// u: array of chars containing boolean values
//static unsigned char *pack_bit_data(unsigned char *u, int nu)
//{
//	int np = nu / 8;
//	if (np*8 != nu)
//		fail("can not unpack an odd number of bits");
//	unsigned char *p = xmalloc(np);
//	for (int i = 0; i < np; i++)
//	{
//		p[i] = 0;
//		for (int j = 0; j < 8; j++)
//			if (u[8*i+j])
//				SETBIT(p[i], j);
//	}
//	return p;
//}
//
//// p: array of chars to be separated into bits
//static unsigned char *unpack_bit_data(unsigned char *p, int np)
//{
//	int nu = 8*np;
//	unsigned char *u = xmalloc(nu);
//	for (int i = 0; i < np; i++)
//		for (int j = 0; j < 8; j++)
//			u[8*i+j] = GETBIT(p[i], j);
//	return u;
//}


// BINARY RUN-LENGTH ENCODING
// (used as a pre-processing for modified huffman coder)
//
// 1st encoded number: status of first bit
// 2nd encoded number: length of first run of equal bits
// 3rd encoded number: length of second run of equal bits
// ...
//

//#define RLE_MAXRUN 10


//// binary run-length encoding
//static int fax_run_length_encoding(unsigned char *out, unsigned char *in, int n)
//{
//	int r = 0;
//	out[r++] = in[0];
//	if (in[0] != 0 && in[0] != 1) fail("bad bad first bit! (%d)", in[0]);
//	out[r] = 1;
//	for (int i = 1; i < n; i++)
//	{
//		if (in[i] != 0 && in[i] != 1) fail("bad bad %d bit! (%d)", i, in[i]);
//		if (out[r] == UINT8_MAX)
//		{
//			out[++r] = 0;
//			out[++r] = 0;
//		}
//		if (in[i] == in[i-1])
//			out[r] += 1;
//		else
//			out[++r] = 1;
//	}
//	return r+1;
//}
//
//// binary run-length decoding
//static int fax_run_length_decoding(unsigned char *out, unsigned char *in, int n)
//{
//	int curr = in[0], r = 0;
//	out[r] = curr;
//	for (int i = 1; i < n; i++)
//	{
//		for (int j = 0; j < in[i]; j++)
//			out[r++] = curr;
//		curr = !curr;
//	}
//	return r;
//}

#include "dataconv.c"

//static struct bitmap_font unpack_font(struct bitmap_font f)
//{
//	f.data = alloc_and_transform_from_RAW_to_BIT(f.data, f.ndata, &f.ndata);
//	f.packing = UNPACKED;
//	return f;
//}
//
//static struct bitmap_font pack_font(struct bitmap_font f)
//{
//	f.data = alloc_and_transform_from_BIT_to_RAW(f.data, f.ndata, &f.ndata);
//	f.packing = PACKED;
//	return f;
//}

//// PACKED => UNPACKED
//static struct bitmap_font unpack_font(struct bitmap_font *fp)
//{
//	struct bitmap_font fu = *fp;
//	int usize = fu.width * fu.height * fu.number_of_glyphs;
//	int psize = usize / 8;
//	fu.data = unpack_bit_data(fp->data, psize);
//	assert(fp->packing == PACKED);
//	fu.packing = UNPACKED;
//	fu.ndata = usize;
//	return fu;
//}
//
//// UNPACKED => PACKED
//static struct bitmap_font pack_font(struct bitmap_font *fu)
//{
//	struct bitmap_font fp = *fu;
//	int usize = fp.width * fp.height * fp.number_of_glyphs;
//	int psize = usize / 8;
//	fp.data = pack_bit_data(fu->data, usize);
//	assert(fu->packing == UNPACKED);
//	fp.packing = PACKED;
//	fp.ndata = psize;
//	return fp;
//}

static char *packing_string(enum font_data_format p)
{
	switch(p) {
#define casepack(s) case s: return #s
	casepack(PACKED);
	casepack(UNPACKED);
	casepack(RLE);
	casepack(RLET);
	casepack(X85);
	casepack(X85RLE);
	casepack(X85RLET);
	casepack(DIFF);
	casepack(XOR);
	casepack(RLEDIFF);
	casepack(RLEXOR);
	casepack(PCX);
	casepack(RLEPCX);
	casepack(RLEXORPCX);
	casepack(PCXX85);
	casepack(RLEPCXX85);
	casepack(RLEXORPCXX85);
#undef casepack
	default: fail("impossible packing style");
	}
}

static enum font_data_format packing_unstring(char *s)
{
#define casepack(x) if (0 == strcmp(s,#x)) return x
	casepack(PACKED);
	casepack(UNPACKED);
	casepack(RLE);
	casepack(RLET);
	casepack(X85);
	casepack(X85RLE);
	casepack(X85RLET);
	casepack(DIFF);
	casepack(XOR);
	casepack(RLEDIFF);
	casepack(RLEXOR);
	casepack(PCX);
	casepack(RLEPCX);
	casepack(RLEXORPCX);
	casepack(PCXX85);
	casepack(RLEPCXX85);
	casepack(RLEXORPCXX85);
#undef casepack
	fail("unrecognized packing \"%s\"", s);
}

//static struct bitmap_font reformat_font(struct bitmap_font *f,
//		enum font_data_format fmt)
//{
//	//fprintf(stderr, "\treformatting call %s => %s\n", packing_string(f->packing), packing_string(fmt));
//	if (fmt == f->packing) return *f;
//	if (fmt == UNPACKED && f->packing == PACKED) return unpack_font(*f);
//	if (fmt == PACKED   && f->packing == UNPACKED) return pack_font(*f);
//	if (fmt == RLE) {
//		struct bitmap_font fu = reformat_font(f, UNPACKED);
//		struct bitmap_font fr = fu;
//		fr.data = xmalloc(2+fu.ndata);
//		fr.ndata = fax_run_length_encoding(fr.data, fu.data, fu.ndata);
//		//fprintf(stderr, "rle %d => %d {%g} (%g)\n", fu.ndata, fr.ndata, fu.ndata*1.0/fr.ndata,fr.ndata*100.0/fu.ndata);
//		fr.packing = RLE;
//		return fr;
//	}
//	if (f->packing == RLE) {
//		unsigned char *urle = xmalloc(129*f->ndata);
//		int nurle = fax_run_length_decoding(urle, f->data, f->ndata);
//		//fprintf(stderr, "urle %d => %d\n", f->ndata, nurle);
//		struct bitmap_font ftmp = *f;
//		ftmp.data = urle;
//		ftmp.ndata = nurle;
//		ftmp.packing = UNPACKED;
//		return reformat_font(&ftmp, fmt);
//	}
//	fail("unimplemented conversion \"%s\"\n", packing_string(fmt));
//}

static struct bitmap_font reformat_font(struct bitmap_font f,
		enum font_data_format fmt)
{
	fprintf(stderr, "reformat_font(%s -> %s)\n",
			packing_string(f.packing), packing_string(fmt));
	uint8_t *(*transform)(uint8_t *t, int, int*) = NULL;

	if (fmt == f.packing) {
		return f;
	} else if (fmt == UNPACKED && f.packing == PACKED) {
		transform = alloc_and_transform_from_RAW_to_BIT;
	} else if (fmt == PACKED   && f.packing == UNPACKED) {
		transform = alloc_and_transform_from_BIT_to_RAW;
	} else if (fmt == RLE) {
		f = reformat_font(f, UNPACKED);
		transform = alloc_and_transform_from_BIT_to_RLE1;
	} else if (f.packing == RLE && fmt == UNPACKED) {
		transform = alloc_and_transform_from_RLE1_to_BIT;
	} else if (f.packing == PCXX85 && fmt == UNPACKED) {
		// PCXX85 -x85toraw-> PCX -pcxtoraw-> PACKED -rawtobit-> UNPACKED
	// TODO arrays of transforms (turn this function from code to data)
	f.data = alloc_and_transform_from_X85_to_RAW(f.data, f.ndata, &f.ndata);
	f.data = alloc_and_transform_from_RLE8_to_RAW(f.data, f.ndata, &f.ndata);
	f.data = alloc_and_transform_from_RAW_to_BIT(f.data, f.ndata, &f.ndata);
	f.packing = UNPACKED;
	return f;
	} else if (fmt == DIFF) {
		f = reformat_font(f, PACKED);
		transform = alloc_and_transform_diff;
	} else if (fmt == XOR) {
		f = reformat_font(f, PACKED);
		transform = alloc_and_transform_xor;
	} else if (fmt == RLEDIFF) {
		f = reformat_font(f, RLE);
		transform = alloc_and_transform_diff;
	} else if (fmt == RLEXOR) {
		f = reformat_font(f, RLE);
		transform = alloc_and_transform_xor;
	} else if (fmt == PCX) {
		f = reformat_font(f, PACKED);
		transform = alloc_and_transform_from_RAW_to_RLE8;
	} else if (fmt == RLEPCX) {
		f = reformat_font(f, RLE);
		transform = alloc_and_transform_from_RAW_to_RLE8;
	} else if (fmt == RLEXORPCX) {
		f = reformat_font(f, RLEXOR);
		transform = alloc_and_transform_from_RAW_to_RLE8;
	} else if (fmt == X85) {
		f = reformat_font(f, PACKED);
		transform = alloc_and_transform_from_RAW_to_X85;
	} else if (fmt == X85RLE) {
		f = reformat_font(f, RLE);
		transform = alloc_and_transform_from_RAW_to_X85;
	} else if (fmt == PCXX85) {
		f = reformat_font(f, PCX);
		transform = alloc_and_transform_from_RAW_to_X85;
	} else if (fmt == RLEPCXX85) {
		f = reformat_font(f, RLEPCX);
		transform = alloc_and_transform_from_RAW_to_X85;
	} else if (fmt == RLEXORPCXX85) {
		f = reformat_font(f, RLEXORPCX);
		transform = alloc_and_transform_from_RAW_to_X85;
	} else if (fmt == RLET) {
		f = reformat_font(f, UNPACKED);
		f.data = alloc_and_transpose_3d1(f.data, f.width, f.height,
				f.number_of_glyphs);
		transform = alloc_and_transform_from_BIT_to_RLE1;
	} else
		fail("unimplemented conversion \"%s\"=>\"%s\"\n",
				packing_string(f.packing), packing_string(fmt));

	if (transform)
		f.data = transform(f.data, f.ndata, &f.ndata);
	f.packing = fmt;
	return f;

	////fprintf(stderr, "\treformatting call %s => %s\n", packing_string(f->packing), packing_string(fmt));
	//if (fmt == f->packing) return *f;
	//if (fmt == UNPACKED && f->packing == PACKED) return unpack_font(*f);
	//if (fmt == PACKED   && f->packing == UNPACKED) return pack_font(*f);
	//if (fmt == RLE) {
	//	struct bitmap_font fu = reformat_font(f, UNPACKED);
	//	struct bitmap_font fr = fu;
	//	fr.data = xmalloc(2+fu.ndata);
	//	fr.ndata = fax_run_length_encoding(fr.data, fu.data, fu.ndata);
	//	//fprintf(stderr, "rle %d => %d {%g} (%g)\n", fu.ndata, fr.ndata, fu.ndata*1.0/fr.ndata,fr.ndata*100.0/fu.ndata);
	//	fr.packing = RLE;
	//	return fr;
	//}
	//if (f->packing == RLE) {
	//	unsigned char *urle = xmalloc(129*f->ndata);
	//	int nurle = fax_run_length_decoding(urle, f->data, f->ndata);
	//	//fprintf(stderr, "urle %d => %d\n", f->ndata, nurle);
	//	struct bitmap_font ftmp = *f;
	//	ftmp.data = urle;
	//	ftmp.ndata = nurle;
	//	ftmp.packing = UNPACKED;
	//	return reformat_font(&ftmp, fmt);
	//}
	fail("unimplemented conversion \"%s\"\n", packing_string(fmt));
}

// this function prints a string that must be a valid rvalue
// for the lvalue "unsigned char x[] = "
static void dump_data_in_c(FILE *f, void *data, int n, enum font_data_format p)
{
	unsigned char *udata = data;
	if (p!=X85 && p!=X85RLE && p!=PCXX85 && p!=RLEPCXX85 && p!=RLEPCXX85)
		fprintf(f, " {\n\t");
	switch(p) {
	case UNPACKED:
		for (int i = 0; i < n; i++) {
			int bit = udata[i];
			if (bit != 0 && bit != 1) fail("bad bit %d\n", bit);
			fprintf(f, bit?"1,%s":"0,%s",(i&&!(i%16))?"\n\t":" ");
		}
		break;
	case PACKED:
		for (int i = 0; i < n; i++)
			fprintf(f, "%d,%s", udata[i], (i&&!(i%16))?"\n\t":" ");
		break;
	case RLE:
	case RLET:
	case DIFF:
	case XOR:
	case RLEDIFF:
	case RLEXOR:
	case RLEPCX:
	case RLEXORPCX:
	case PCX:
		for (int i = 0; i < n; i++)
			fprintf(f, "%d,%s", udata[i], (i&&!(i%16))?"\n\t":" ");
		break;
	case X85:
	case X85RLE:
	case PCXX85:
	case RLEPCXX85:
	case RLEXORPCXX85:
		{
			fprintf(f, "\n");
			int llen = 77, ns = n/llen;
			for (int i = 0; i < ns; i++)
			{
				fprintf(f, "\"");
				for (int j = 0; j < llen; j++)
					fprintf(f, "%c", udata[i*llen+j]);
				fprintf(f, "\"\n");
			}
			fprintf(f, "\"");
			for (int j = ns*llen; j < n; j++)
				fprintf(f, "%c", udata[j]);
			fprintf(f, "\"");

		}
		break;
	default: fail("bad dump format by now");
	}
	if (p!=X85 && p!=X85RLE && p!=PCXX85 && p!=RLEPCXX85 && p!=RLEPCXX85)
		fprintf(f, "0\n}");
}

static void dump_font_as_parseable_c_struct_as_it_is(char *filename,
		struct bitmap_font font, char *name)
{
	FILE *f = xfopen(filename, "w");
	fprintf(f, "unsigned char %s_data[] =", name);
	dump_data_in_c(f, font.data, font.ndata, font.packing);
	fprintf(f, ";\n");
	fprintf(f, "struct bitmap_font %s[1] = "
			"{{%d, %d, %d, %s, \"%s\", %s_data}};\n",
		name, font.number_of_glyphs, font.width, font.height,
		packing_string(font.packing), name, name);
	xfclose(f);

	double e = entropy(font.data, font.ndata);
	double cs = font.ndata * e/8;
	double cf = 100.0*cs/font.ndata;
	fprintf(stderr, "entropy %s %s  = %g  %d => %g %g%%\n",
			name, packing_string(font.packing),
			e,font.ndata, cs, cf);
}

static void dump_font_as_parseable_c_struct(char *filename,
		struct bitmap_font font, char *name)
{
	struct bitmap_font cfont = reformat_font(font, PACKED);
	dump_font_as_parseable_c_struct_as_it_is(filename, cfont, name);
}

//static void dump_font_as_parseable_c_struct_old(char *filename,
//		struct bitmap_font *ufont, char *name)
//{
//	struct bitmap_font font = pack_font(ufont);
//	FILE *f = xfopen(filename, "w");
//	int datasize = font.number_of_glyphs * font.width * font.height;
//	fprintf(f, "unsigned char %s_data[] = {\n\t", name);
//	for (int i = 0; i < datasize/8; i++)
//		fprintf(f, "%d,%s", (unsigned char)font.data[i],
//				(i&&!(i%16))?"\n\t":" ");
//	fprintf(f, "0\n};\n");
//	fprintf(f, "struct bitmap_font %s[1] = "
//			"{{%d, %d, %d, PACKED, %s_data}};\n",
//		name, font.number_of_glyphs, font.width, font.height, name);
//	xfclose(f);
//	free(font.data);
//}


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
	//assert(i >= 0);
	//assert(j >= 0);
	//assert(i < f->width);
	//assert(j < f->height);
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
	return sscanf(line, " %d %d %d %d\n", nums, nums+1, nums+2, nums+3);
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
		if (r != 1) fail("could not read BDF tag STARTFONT\n");
		//fprintf(stderr, "BDF version \"%g\"\n", bdfversion);
	}
	font->data = NULL;
	font->number_of_glyphs = NGLYPHS;
	font->width = font->height = -1;
	int bitmap_line = -2, glyph = -1, offx=0, offy=0, descent=0;
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
			offx = 0;//num[2];
			offy = 0;//num[3];
			size_t ds = num[0] * num[1] * NGLYPHS;
			font->data = xmalloc(ds);
			font->packing = UNPACKED;
			font->ndata = ds;
			memset(font->data, 0, ds);
		}
		//if (0 == strcmp(tag, "FONT_ASCENT")) ascent = num[0];
		if (0 == strcmp(tag, "FONT_DESCENT")) descent = num[0];
		//if (0 == strcmp(tag, "STARTCHAR")) fprintf(stderr, "BDF \"%s\"\n", line);
		if (0 == strcmp(tag, "ENCODING") && nn == 1) glyph = num[0];
		if (0 == strcmp(tag, "BITMAP")) bitmap_line = -1;
		if (0 == strcmp(tag, "ENDCHAR")) bitmap_line = -2;
		if (0 == strcmp(tag, "BBX") && nn == 4)
		{
			offx = num[2];//font->width - num[0] - num[2];
			offy = font->height - num[1] - num[3] - descent;
			//if (islower(glyph)){
			//fprintf(stderr, "bbx tag %d %d %d %d\t", num[0], num[1], num[2], num[3]);
			//fprintf(stderr, "off glyph %d '%c' = %d %d\n", glyph, glyph, offx, offy);
			//}
		}
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
				{
					int ii = i + offx;
					int jj = bitmap_line + offy;
					set_font_bit(font, glyph, ii, jj);
				}
		}
		if (bitmap_line > -2) bitmap_line += 1;
	}
	xfclose(f);
	if (!font->data)
		fail("failed to parse BDF file \"%s\"\n", fname);
}

static void put_pixel(float *x, int w, int h, int pd, int i, int j, float *c)
{
	//if (j*w + i < w*h)
	if (i>=0 && j>=0 && i<w && j<h)
		for (int l = 0; l < pd; l++)
			x[(w*j+i)*pd+l] = c[l];
}

static void put_pixel_rgb(uint8_t *x, int w, int h, int i, int j, uint8_t *c)
{
	if (c && i>=0 && j>=0 && i<w && j<h)
		for (int l = 0; l < 3; l++)
			x[(w*j+i)*3+l] = c[l];
}

static void put_string_in_float_image(float *x, int w, int h, int pd,
		int posx, int posy, float *color, int kerning,
		struct bitmap_font *font, char *string)
{
	int posx0 = posx;
	while (1)
	{
		int c = *string++;
		if (!c) break;
		if (c == '\n') {
			posx = posx0;
			posy += font->height;
			continue;
		}
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
			       	//else {
				//	float white[10] = {255};
				//	int ii = posx + i;
				//	int jj = posy + j;
				//	put_pixel(x, w, h, pd, ii, jj, white);
				//}
		}
		posx += font->width + kerning;
	}
}

static void put_string_in_rgb_image(uint8_t *x, int w, int h,
		int posx, int posy, uint8_t *fg, uint8_t *bg, int kerning,
		struct bitmap_font *font, char *string)
{
	int posx0 = posx;
	while (1)
	{
		int c = *string++;
		if (!c) break;
		if (c == '\n') {
			posx = posx0;
			posy += font->height;
			continue;
		}
		if (c > 0 && c < font->number_of_glyphs)
		{
			for (int i = 0; i < font->width; i++)
			for (int j = 0; j < font->height; j++)
			{
				int ii = posx + i;
				int jj = posy + j;
				if (get_font_bit(font, c, i, j))
					put_pixel_rgb(x, w, h, ii, jj, fg);
				else
					put_pixel_rgb(x, w, h, ii, jj, bg);
			}
		}
		posx += font->width + kerning;
	}
}

#ifndef OMIT_MAIN_FONTU
#define MAIN_FONTU
#endif//OMIT_MAIN_FONTU

#ifdef MAIN_FONTU

// fontu cdump name {packed|unpacked|zrle|huffman|lzw|mhuffman} [in.bdf [out.c]]
// fontu puts [-f bdf] [-c color] "string" [in.png [out.png]]

static int main_cdump(int c, char **v)
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t"
			"%s name [in.bdf [out.c]]\n", *v);
		//        0  1    2       3
		return 1;
	}
	char *name = v[1];
	char *filename_in = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	struct bitmap_font f;
	font_fill_from_bdf(&f, filename_in);
	dump_font_as_parseable_c_struct(filename_out, f, name);

	return 0;
}

static int main_cdumpf(int c, char **v)
{
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
			"%s name {packing} [in.bdf [out.c]]\n", *v);
		//        0  1    2         3       4
		return 1;
	}
	char *name = v[1];
	enum font_data_format fmt = packing_unstring(v[2]);
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	struct bitmap_font f;
	font_fill_from_bdf(&f, filename_in);
	f = reformat_font(f, fmt);
	dump_font_as_parseable_c_struct_as_it_is(filename_out, f, name);

	free(f.data);

	return 0;
}


static int main_dumptry(int c, char **v)
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s in.bdf\n", *v);
		return 1;
	}
	struct bitmap_font f;
	font_fill_from_bdf(&f, v[1]);

	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout1.c", f, "nom");

	fprintf(stderr, "formatting to packed\n");
	f = reformat_font(f, PACKED);
	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout2.c", f, "nom");

	fprintf(stderr, "formatting to unpacked\n");
	f = reformat_font(f, UNPACKED);
	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout3.c", f, "nom");

	fprintf(stderr, "formatting to packed\n");
	f = reformat_font(f, PACKED);
	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout4.c", f, "nom");

	fprintf(stderr, "formatting to rle\n");
	f = reformat_font(f, RLE);
	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout5.c", f, "nom");

	fprintf(stderr, "formatting to packed\n");
	f = reformat_font(f, PACKED);
	dump_font_as_parseable_c_struct_as_it_is("/tmp/fontout6.c", f, "nom");

	return 0;
}

// fontu puts [-f bdf] [-c color] "string" [in.png [out.png]]

#include "iio.h"

#include "pickopt.c"

static int main_puts(int c, char **v)
{
	int kerning = atoi(pick_option(&c, &v, "k", "0"));
	char *colorname = pick_option(&c, &v, "c", "000");
	if (c != 5 && c != 6 && c != 7) {
		fprintf(stderr, "usage:\n\t"
			"%s font.bdf px py \"string\" [in [out]]\n", *v);
		//        0  1       2  3    4         5   6
		return 1;
	}
	char *bdf = v[1];
	int px = atoi(v[2]);
	int py = atoi(v[3]);
	char *text = v[4];
	char *filename_in = c > 5 ? v[5] : "-";
	char *filename_out = c > 6 ? v[6] : "-";

	struct bitmap_font f;
	font_fill_from_bdf(&f, bdf);

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	float color[10] = {0};
	if (pd == (int)strlen(colorname))
		for (int i = 0; i < pd; i++)
			color[i] = (unsigned char)((255*(colorname[i]-'0'))/8);
	put_string_in_float_image(x, w,h,pd, px,py, color, kerning, &f, text);

	iio_save_image_float_vec(filename_out, x, w, h, pd);

	free(f.data);
	free(x);
	return 0;
}

int main(int c, char **v)
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "cdump")) return main_cdump(c-1, v+1);
	else if (0 == strcmp(v[1], "cdumpf")) return main_cdumpf(c-1, v+1);
	else if (0 == strcmp(v[1], "dumptry")) return main_dumptry(c-1, v+1);
	else if (0 == strcmp(v[1], "puts")) return main_puts(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s {cdump|puts} params\n", *v);
	       return 1;
	}
}
#endif//MAIN_FONTU
