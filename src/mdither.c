// mdither
//
// program for encoding/decoding data into binary dithering patterns
//
//
// mdither count in.png                   # print the image bit capacity
// mdither encode in.png out.png < bytes  # prints the number of bits encoded
// mdither decode in.png > bytes




// 0   8  4  2  1   12 10 9  6  5  3   7  11 13 14  15
//
// 00  10 01 00 00  11 10 10 01 01 00  01 10 11 11  11
// 00  00 00 10 01  00 10 01 10 01 11  11 11 01 10  11
//
// X   a  a  b  b   X  c  c  d  d  X   e  e  f  f   X

// what bit is encoded by each pattern
const int table_decoder_2x2[16] = {
	[0]  = -1, [15] = -1,
	[12] = -1, [3]  = -1,
	[8]  = 0,  [4]  = 1,
	[2]  = 0,  [1]  = 1,
	[10] = 0,  [9]  = 1,
	[6]  = 0,  [5]  = 1,
	[7]  = 0,  [11] = 1,
	[13] = 0,  [14] = 1,
};

// inverse of the previous table
const int table_encoder_2x2[16][2] = {
	[0]  = {-1, -1},
	[1]  = { 2,  1},
	[2]  = { 2,  1},
	[3]  = {-1, -1},
	[4]  = { 8,  4},
	[5]  = { 6,  5},
	[6]  = { 6,  5},
	[7]  = { 7, 11},
	[8]  = { 8,  4},
	[9]  = {10,  9},
	[10] = {10,  9},
	[11] = { 7, 11},
	[12] = {-1, -1},
	[13] = {13, 14},
	[14] = {13, 14},
	[15] = {-1, -1},
};

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>


// extract the bits encoded in a pattern
// return value: number of decoded bits
int cell_decode(
		int *o,  // table of decoded bits (filled-in by this algorithm)
		int p    // input bit pattern (indexed)
		)
{
	assert(p >= 0);
	assert(p <= 15);
	int b = table_decoder_2x2[p];
	if (b < 0)
		return 0;
	if (o)
		o[0] = b;
	return 1;
}

// encode a bit into a pattern
// return value: the encoded pattern, -1 if encoding is not possible
int cell_encode_one_bit(
		int p, // a pattern before encoding
		int b  // the bit to encode
		)
{
	assert(p >= 0 && p <= 15);
	assert(b==0 || b==1);
	return table_encoder_2x2[p][b];
}

// decode all the cells from a binary image
// return value: number of decoded bits
int decode_image(
		int *B,  // table of bits to be filled-in
		float *x, // input image data (boolean given by x>0)
		int w,    // input image width
		int h     // input image width
		)
{
	int r = 0; // index in the output bit table
	for (int j = 0; j < h-1; j += 2)
	for (int i = 0; i < w-1; i += 2)
	{
		// a b    8 4
		// c d    2 1
		int a = x[(j+0)*w+(i+0)] > 0;
		int b = x[(j+0)*w+(i+1)] > 0;
		int c = x[(j+1)*w+(i+0)] > 0;
		int d = x[(j+1)*w+(i+1)] > 0;
		int p = d + 2*c + 4*b + 8*a;
		r += cell_decode(B + r, p);
	}
	return r;
}

// encode zeros and ones from stdin into the given image
// return value: number of lines read from stdin
int encode_bits_into_image(
		float *y,     // output image data
		float *x,     // input image data
		int w,        // image width
		int h,        // image height
		int *B,       // input bits (ints equal to 0 or one)
		int n         // maximum number of input bits
		)
{
	for (int i = 0; i < w*h; i++)
		y[i] = x[i];

	int r = 0; // number of bits encoded
	for (int j = 0; j < h-1; j += 2)
	for (int i = 0; i < w-1; i += 2)
	{
		// a b    8 4
		// c d    2 1
		int a = x[(j+0)*w+(i+0)] > 0;
		int b = x[(j+0)*w+(i+1)] > 0;
		int c = x[(j+1)*w+(i+0)] > 0;
		int d = x[(j+1)*w+(i+1)] > 0;
		int p = d + 2*c + 4*b + 8*a;
		assert(p >= 0 && p <= 15);
		assert(B[r] == 0 || B[r] == 1);
		int q = table_encoder_2x2[p][B[r]];
		if (q < 0 || r+1 >= n) { // no encoding possible, just copy p
			y[(j+0)*w+(i+0)] = a;
			y[(j+0)*w+(i+1)] = b;
			y[(j+1)*w+(i+0)] = c;
			y[(j+1)*w+(i+1)] = d;
		} else { // encode the bits of q
			y[(j+0)*w+(i+0)] = (q/8)%2;
			y[(j+0)*w+(i+1)] = (q/4)%2;
			y[(j+1)*w+(i+0)] = (q/2)%2;
			y[(j+1)*w+(i+1)] = (q/1)%2;
			r += 1;
		}
	}
	return r;
}

// query how many bits can be encoded in this particular binary image
int image_capacity_in_bits(
		float *x, // input image data (boolean given by x>0)
		int w,    // input image width
		int h     // input image width
		)
{
	int r = 0; // bit counter
	for (int j = 0; j < h-1; j += 2)
	for (int i = 0; i < w-1; i += 2)
	{
		// a b    8 4
		// c d    2 1
		int a = x[(j+0)*w+(i+0)] > 0;
		int b = x[(j+0)*w+(i+1)] > 0;
		int c = x[(j+1)*w+(i+0)] > 0;
		int d = x[(j+1)*w+(i+1)] > 0;
		int p = d + 2*c + 4*b + 8*a;
		r += cell_decode(NULL, p);
	}
	return r;
}



// auxiliary function to write bits into a file
static void write_n_bits_to_file(FILE *f, int *b, int n)
{
	int c = 0;
	for (int i = 0; i < n; i++)
	{
		c += b[i] << (i%8);
		if (7 == i%8 || i==n-1)
		{
			fputc(c, f);
			c = 0;
		}
	}
}

// auxiliary function to read bits from a file
static void read_n_bits_from_file(int *b, int n, FILE *f)
{
	for (int i = 0; i < n; i++)
		b[i] = 0;
	int c = 0;
	for (int i = 0; i < n; i++)
	{
		if (0 == i%8)
			c = fgetc(f);
		if (c == EOF)
			break;
		b[i] = (bool) (c & ( 1 << (i%8) ) );
	}
}


#include "iio.h"


// CLI for debugging purposes (pipes 3 bytes from stdin to stdout)
int main_bcat(int c, char *v[])
{
	int b[24];
	read_n_bits_from_file(b, 24, stdin);
	for (int i = 0; i < 24; i++)
		fprintf(stderr, "b[%d] = %d\n", i, b[i]);
	write_n_bits_to_file(stdout, b, 24);
	return 0;
}


// CLI for decoding
int main_mddecode(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s img > bits.txt\n", *v);
		//                                 0 1
	char *filename_in  = v[1];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	int *b = malloc(w*h);
	int n = decode_image(b, x, w, h);
	fprintf(stderr, "decoded %d bits from an image of sixe %dx%d\n", n,w,h);
	write_n_bits_to_file(stdout, b, n);
	return 0;
}


// CLI for encoding
int main_mdencode(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s in out < bytes\n", *v);
		//                                 0 1  2
	char *filename_in  = v[1];
	char *filename_out = v[2];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	float *y = malloc(w*h*sizeof*x);
	int n = image_capacity_in_bits(x, w, h);
	fprintf(stderr, "cap(%d,%d)=%d\n", w, h, n);
	int *b = malloc(n*sizeof*b);
	read_n_bits_from_file(b, n, stdin);
	int m = encode_bits_into_image(y, x, w, h, b, n);
	fprintf(stderr, "encoded %d bits into an image of size %dx%d\n",m,w,h);
	for (int i = 0; i < w*h; i++)
		y[i] *= 255;
	iio_write_image_float(filename_out, y, w, h);
	return 0;
}


// CLI for counting the bit capacity
int main_mdcount(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s image.png\n", *v);
	//                                         0 1
	char *filename_in = v[1];
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	int n = image_capacity_in_bits(x, w, h);
	printf("%d bits %g bytes %g k %g M\n", n, n/8.0, n/8192.0, n/8388608.0);
	return 0;
}

// main selector
int main_mdither(int c, char *v[])
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "count")) return main_mdcount(c-1, v+1);
	else if (0 == strcmp(v[1], "encode")) return main_mdencode(c-1, v+1);
	else if (0 == strcmp(v[1], "decode")) return main_mddecode(c-1, v+1);
	else if (0 == strcmp(v[1], "bcat")) return main_bcat(c-1, v+1);
	else {
	usage:
		fprintf(stderr, "usage:\n\t%s [count|encode|decode] ...\n", *v);
		return 0;
	}
}


#ifndef HIDE_ALL_MAINS
int main(int c, char **v){ return main_mdither(c, v); }
#endif
