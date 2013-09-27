// BINARY DATA FORMATS
// -------------------
//
// RAW: bytes, or "packed bits"
// BIT: one bit on every byte, or "unpacked bits"
// RLE1: bits encoded into bytes by lenght of constant runs
// RLE8: run-length encoded bytes, à la PCX
// HUF8: canonically Huffman-encoded bytes
// HUF1: canonically Huffman-encoded runs of bits (i.e., RLE1 followed by HUF8)
// B64: base 64 encoding of bytes into printable characters
// A85: base 85 encoding of bytes into printable characters
// X85: base 85 encoding of bytes into printable characters (alternate charset)
// FAX: RLE1 + HUF8 + B64
// LZW: Lempel-Ziv (not implemented)

// implementation idea: do not think in terms of data conversions between pairs
// of formats, but in terms of data transformations, from whatever format into
// the target one.  That way, composition of formats is automatically supported.


///////////////
////  API  ////
///////////////
#include <stdint.h>


// 1. elementary transformations

// 1.1. unpacking and packing
uint8_t *alloc_and_transform_from_RAW_to_BIT(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_BIT_to_RAW(uint8_t *x, int n, int *nout);

// 1.2. run-length encoding and decoding
uint8_t *alloc_and_transform_from_BIT_to_RLE1(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_RLE1_to_BIT(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_RAW_to_RLE8(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_RLE8_to_RAW(uint8_t *x, int n, int *nout);

// 1.3. differential coding and decoding
uint8_t *alloc_and_transform_diff(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_undiff(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_xor(uint8_t *x, int n, int *nout);

// 1.3. canonical huffman encoding and decoding
uint8_t *alloc_and_transform_from_RAW_to_HUF8(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_HUF8_to_RAW(uint8_t *x, int n, int *nout);

// 1.4. base64 encoding
uint8_t *alloc_and_transform_from_RAW_to_B64(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_B64_to_RAW(uint8_t *x, int n, int *nout);

// 1.5. base85 encoding
uint8_t *alloc_and_transform_from_RAW_to_A85(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_A85_to_RAW(uint8_t *x, int n, int *nout);

// 1.6. lempel-ziv compression (will not implemented)
//uint8_t *alloc_and_transform_from_RAW_to_LZW(uint8_t *x, int n, int *nout);



// 2. composite transformations (as a form of API convenience)

// 2.1. BIT => RLE1 => HUF8
uint8_t *alloc_and_transform_from_BIT_to_HUF1(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_HUF1_to_BIT(uint8_t *x, int n, int *nout);

// 2.2. RAW => BIT => RLE1 => HUF8
uint8_t *alloc_and_transform_from_RAW_to_HUF1(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_HUF1_to_RAW(uint8_t *x, int n, int *nout);

// 2.3. BIT => RLE1 => HUF8 => B64
uint8_t *alloc_and_transform_from_BIT_to_B64(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_B64_to_BIT(uint8_t *x, int n, int *nout);

// 2.4. RAW => BIT => RLE1 => HUF8 => B64
uint8_t *alloc_and_transform_from_RAW_to_FAX(uint8_t *x, int n, int *nout);
uint8_t *alloc_and_transform_from_FAX_to_RAW(uint8_t *x, int n, int *nout);



// 3. ordering transforms
uint8_t *alloc_and_transpose_2d(uint8_t *x, int w, int h);
uint8_t *alloc_and_transpose_3dx(uint8_t *x, int w, int h, int d);
uint8_t *alloc_and_transpose_3dy(uint8_t *x, int w, int h, int d);
uint8_t *alloc_and_transpose_3dz(uint8_t *x, int w, int h, int d);


//////////////////////////
////  IMPLEMENTATION  ////
//////////////////////////

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "xmalloc.c"
#include "fail.c"

#define SETBIT(x,i) ((x)|=(1<<(i)))
#define GETBIT(x,i) (bool)((x)&(1<<(i)))


// unpack bytes into individual bits, thus enlarging the array eightfold
uint8_t *alloc_and_transform_from_RAW_to_BIT(uint8_t *x, int n, int *nout)
{
	*nout = 8*n;
	uint8_t *y = xmalloc(*nout+16);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < 8; j++)
			y[8*i+j] = GETBIT(x[i], j);
	return y;
}

// pack groups of 8 bits into bytes, thus reducing the array size eightfold
uint8_t *alloc_and_transform_from_BIT_to_RAW(uint8_t *x, int n, int *nout)
{
	*nout = n / 8;
	if (*nout * 8 != n)
		fail("can not unpack an odd number (%d) of bits", n);
	uint8_t *y = xmalloc(*nout+1);
	for (int i = 0; i < *nout; i++)
	{
		y[i] = 0;
		for (int j = 0; j < 8; j++)
			if (x[8*i+j])
				SETBIT(y[i], j);
	}
	return y;
}

// run-length encode a sequence of bits into runs
uint8_t *alloc_and_transform_from_BIT_to_RLE1(uint8_t *x, int n, int *nout)
{
	int r = 0;
	uint8_t *y = xmalloc(n+25);
	y[r++] = (bool)x[0];
	y[r] = 1;
	for (int i = 1; i < n; i++)
	{
		if (y[r] == UINT8_MAX)
		{
			y[++r] = 0;
			y[++r] = 0;
		}
		if ((bool)x[i] == (bool)x[i-1])
			y[r] += 1;
		else
			y[++r] = 1;
	}
	*nout = r + 1;
	for (int i = *nout; i < *nout+24; i++) y[i] = 0;
	return y;
}

// run-length decoding
uint8_t *alloc_and_transform_from_RLE1_to_BIT(uint8_t *x, int n, int *nout)
{
	bool curr = x[0];
	int r = 0;
	uint8_t *y = xmalloc(n*UINT8_MAX); // worst case: block-constant seq
	y[r] = curr;
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < x[i]; j++)
			y[r++] = curr;
		curr = !curr;
	}
	*nout = r;
	return y;
}

static uint8_t encode_b64_quark(uint8_t x)
{
	assert(x < 64);
	if (x < 26) return 'A' + x;
	if (x < 52) return 'a' + x - 26;
	if (x < 62) return '0' + x - 52;
	if (x == 62) return '+';
	if (x == 63) return '/';
	fail("impossible base64 quark");
}

static uint8_t decode_b64_quark(uint8_t x)
{
	if (isupper(x)) return x - 'A';
	if (islower(x)) return x - 'a' + 26;
	if (isdigit(x)) return x - '0' + 52;
	if (x == '+' || x == '-') return 62;
	if (x == '/' || x == '_') return 63;
	fail("bad quark %c", x);
}

// base 64 encoding
uint8_t *alloc_and_transform_from_RAW_to_B64(uint8_t *x, int n, int *nout)
{
	int nbits;
	uint8_t *bits = alloc_and_transform_from_RAW_to_BIT(x, n, &nbits);
	assert(0 == nbits % 8);
	// TODO: replace these conditions by a single arithmetic expression
	int nrem;
	{	if (0 == nbits % 6) nrem = 0;
		else if (0 == (nbits+8) % 6) nrem = 8;
		else if (0 == (nbits+16) % 6) nrem = 16;
		else fail("impossible condition nbits=%d", nbits);
	}
	for (int i = nbits; i < nbits+nrem; i++)
		bits[i] = 0;
	nbits += nrem;
	assert(0 == nbits % 6);
	int nsix = nbits/6;
	uint8_t *y = xmalloc(nsix);
	for (int i = 0; i < nsix; i++)
	{
		uint8_t six = 0;
		for (int j = 0; j < 6; j++)
			if (bits[6*i+j])
				SETBIT(six, j);
		y[i] = encode_b64_quark(six);
	}
	free(bits);
	*nout = nsix;
	return y;
}

// base 64 decoding
uint8_t *alloc_and_transform_from_B64_to_RAW(uint8_t *x, int n, int *nout)
{
	uint8_t *tmp = xmalloc(n);
	for (int i = 0; i < n; i++)
		tmp[i] = decode_b64_quark(x[i]);

	int nbits;
	uint8_t *bits = alloc_and_transform_from_RAW_to_BIT(tmp, n, &nbits);

	int ibit = 0;
	for (int i = 0; i < nbits; i++)
		if (i%8 < 6)
			bits[ibit++] = bits[i];

	assert(4*ibit == 3*nbits);

	uint8_t *y = alloc_and_transform_from_BIT_to_RAW(bits, ibit, nout);
	free(tmp);
	free(bits);
	return y;
}

// base 85 encoding
uint8_t *alloc_and_transform_from_RAW_to_A85(uint8_t *x, int n, int *nout)
{
	assert(0 == n % 4);
	*nout = 5*n/4;
	uint8_t *y = xmalloc(*nout);
	for (int i = 0; i < n/4; i++)
	{
		uint32_t w = ((uint32_t*)x)[i];
		for (int j = 0; j < 5; j++)
		{
			y[5*i+4-j] = '!' + w % 85;
			w /= 85;
		}
	}
	return y;
}

// base 85 decoding
uint8_t *alloc_and_transform_from_A85_to_RAW(uint8_t *x, int n, int *nout)
{
	assert(0 == n % 5);
	*nout = 4*n/5;
	uint32_t *y = xmalloc(*nout);
	for (int i = 0; i < n/5; i++)
	{
		y[i] = 0;
		for (int j = 0; j < 5; j++)
			y[i] = y[i]*85 + (x[5*i+j] - '!');
	}
	return (void*)y;
}

// base 85 encoding with quotable charset
uint8_t *alloc_and_transform_from_RAW_to_X85(uint8_t *x, int n, int *nout)
{
	if (0 != n % 4)
		n += 4 - n % 4;
	assert(0 == n % 4);
	*nout = 5*n/4;
	uint8_t *y = xmalloc(*nout);
	for (int i = 0; i < n/4; i++)
	{
		uint32_t w = ((uint32_t*)x)[i];
		for (int j = 0; j < 5; j++)
		{
			uint8_t q = w % 85;
			y[5*i+4-j] = q + (q < 52 ? 40 : 41);
			w /= 85;
		}
	}
	return y;
}

// base 85 decoding with quotable charset
uint8_t *alloc_and_transform_from_X85_to_RAW(uint8_t *x, int n, int *nout)
{
	assert(0 == n % 5);
	*nout = 4*n/5;
	uint32_t *y = xmalloc(*nout);
	for (int i = 0; i < n/5; i++)
	{
		y[i] = 0;
		for (int j = 0; j < 5; j++)
		{
			uint8_t q = x[5*i+j];
			assert(isprint(q));
			y[i] = y[i]*85 + q - (q < 92 ? 40 : 41);
		}
	}
	return (void*)y;
}

uint8_t *alloc_and_transpose_2d(uint8_t *x, int w, int h)
{
	uint8_t *y = xmalloc(w*h);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
		y[i*h+j] = x[j*w+i];
	return y;
}

uint8_t *alloc_and_transpose_3d3(uint8_t *x, int w, int h, int d)
{
	uint8_t *y = xmalloc(w*h*d);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	for (int k = 0; k < d; k++)
		// swap the roles of (i,j)
		y[(i*h+j)*d+k] = x[(j*w+i)*d+k];
	return y;
}

uint8_t *alloc_and_transpose_3d2(uint8_t *x, int w, int h, int d)
{
	uint8_t *y = xmalloc(w*h*d);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	for (int k = 0; k < d; k++)
		// swap the roles of (i,k)
		y[(j*d+k)*w+i] = x[(j*w+i)*d+k];
	return y;
}

uint8_t *alloc_and_transpose_3d1(uint8_t *x, int w, int h, int d)
{
	uint8_t *y = xmalloc(w*h*d);
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
	for (int k = 0; k < d; k++)
		// swap the roles of (j,k)
		y[(k*w+i)*h+j] = x[(j*w+i)*d+k];
	return y;
}

uint8_t *alloc_and_transform_diff(uint8_t *x, int n, int *nout)
{
	uint8_t *y = xmalloc(n);
	*nout = n;
	*y = *x;
	for (int i = 1; i < n; i++)
		y[i] = x[i] - x[i-1];
	return y;
}

uint8_t *alloc_and_transform_undiff(uint8_t *x, int n, int *nout)
{
	uint8_t *y = xmalloc(n);
	*nout = n;
	*y = *x;
	for (int i = 1; i < n; i++)
		y[i] = x[i] + x[i-1];
	return y;
}

uint8_t *alloc_and_transform_xor(uint8_t *x, int n, int *nout)
{
	uint8_t *y = xmalloc(n);
	*nout = n;
	*y = *x;
	for (int i = 1; i < n; i++)
		y[i] = x[i] ^ x[i-1];
	return y;
}

// PCX encoding
uint8_t *alloc_and_transform_from_RAW_to_RLE8(uint8_t *x, int n, int *nout)
{
	uint8_t curr = x[0]+1, *y = xmalloc(2*n);
	int r = 0, runlen = 1;
	if (curr == x[1]) curr++;
	assert(x[0] != curr);
	assert(x[1] != curr);
	for (int i = 0; i < n; i++)
	{
		if (i+1 < n && curr == x[i+1] && runlen < 64) // inside run
			runlen += 1;
		else // end of run
		{
			if (x[i] > 191 || runlen > 1)
				y[r++] = 192 + runlen - 1;
			y[r++] = x[i];
			runlen = 1;
			curr = x[i+1];
		}
	}
	*nout = r;
	return y;
}

// PCX decoding
uint8_t *alloc_and_transform_from_RLE8_to_RAW(uint8_t *x, int n, int *nout)
{
	uint8_t *y = xmalloc(64*n);
	int r = 0, i = 0;
	while (i < n)
	{
		if (x[i] > 191) {
			assert(i+1 < n);
			int count = x[i] - 192 + 1;
			for (int j = 0; j < count; j++)
				y[r++] = x[i+1];
			i += 2;
		} else
			y[r++] = x[i++];
	}
	*nout = r;
	return y;
}

double entropy(uint8_t *x, int n)
{
	double t[256] = {0};
	for (int i = 0; i < n; i++)
		t[x[i]] += 1;
	double e = 0;
	for (int i = 0; i < 256; i++)
	{
		double p = t[i] / n;
		//fprintf(stderr, "t[%d] = %g, p=%g\n", i, t[i], p);
		if (p > 0)
			e -= p * log2(p);
	}
	return e;
}



#ifdef MAIN_DATACONV

static void f(void)
{
	//int n = 22;
	//uint8_t data[22] = {1, 2, 3,  7, 7, 7,  1, 2, 3,
	//	200, 2, 3,  222, 222, 222,  1, 2, 3, 44, 45, 46, 47};
	int n = 100003;
	uint8_t data[n];
	for (int i = 0; i < n; i++)
		data[i] = rand();



	fprintf(stderr, "\n\n\noriginal data:\n");
	for (int i = 0; i < n; i++)
		fprintf(stderr, "%d%s", data[i], (i&&(0==(i+1)%77))?"\n":" ");
	fprintf(stderr, "\n");
	fprintf(stderr, "entropy = %g\n", entropy(data, n));


	int ne;
	uint8_t *edata = alloc_and_transform_from_RAW_to_RLE8(data, n, &ne);

	//int ne2;
	//uint8_t *edata2 = alloc_and_transform_from_RAW_to_B64(data, n, &ne2);


	fprintf(stderr, "\nPCX data:\n");
	for (int i = 0; i < ne; i++)
		fprintf(stderr, "%d%s", edata[i], (i&&(0==(i+1)%77))?"\n":" ");
	fprintf(stderr, "\n");
	fprintf(stderr, "%d => %d\n", n, ne);
	fprintf(stderr, "entropy = %g\n", entropy(edata, ne));

	int dne;
	uint8_t *dedata = alloc_and_transform_from_RLE8_to_RAW(edata, ne, &dne);



	fprintf(stderr, "\ndata recovered from PCX:\n");
	for (int i = 0; i < dne; i++)
		fprintf(stderr, "%d%s", dedata[i], (i&&(0==(i+1)%77))?"\n":" ");
	fprintf(stderr, "\n");

	fprintf(stderr, "n=%d dne=%d\n", n, dne);
	assert(n == dne);
	for (int i = 0; i < n; i++)
		assert(data[i] == dedata[i]);


	free(edata);
	free(dedata);
	//free(edata2);
}

int main(void)
{
#define n 24
	uint8_t raw[n] = {
		1, 2, 3, 4, 72, 19, 200, 201,
		202, 0, 9, 100, 120, 130, 140, 127,
		13, 47, 0, 255, 72, 19, 200, 201,
	};

	printf("beginning %d things:\n", n);
	for (int i = 0; i < n; i++)
		printf("%d%c", raw[i], (i&&(0==(i+1)%8))?'\n':' ');

	int nbits;
	uint8_t *bits = alloc_and_transform_from_RAW_to_BIT(raw, n, &nbits);

	//printf("got %d bits:\n", nbits);
	//for (int i = 0; i < nbits; i++)
	//	printf("%d%c", bits[i], (i&&(0==(i+1)%8))?'\n':' ');

	int nb64;
	uint8_t *b64 = alloc_and_transform_from_RAW_to_X85(raw, n, &nb64);
	printf("got %d things:\n", nb64);
	for (int i = 0; i < nb64; i++)
		printf("%c%s", b64[i], (i&&(0==(i+1)%5))?"\n":"");

	int nfin;
	uint8_t *fin = alloc_and_transform_from_X85_to_RAW(b64, nb64, &nfin);
	printf("\ngot %d final things:\n", nfin);
	for (int i = 0; i < nfin; i++)
		printf("%d%c", fin[i], (i&&(0==(i+1)%8))?'\n':' ');

	free(bits);
	free(b64);
	free(fin);

	f();

	return 0;
}


#endif//MAIN_DATACONV
