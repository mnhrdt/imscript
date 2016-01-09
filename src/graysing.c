#include <assert.h>
#include <stdio.h>


//#define SING_REGULAR 123
//#define SING_CYCLIC 132
//#define SING_SINGULAR 231
//#define SING_DEGENERATE_0012 12
//#define SING_DEGENERATE_0021 21
//#define SING_DEGENERATE_0112 112
//#define SING_DEGENERATE_0120 120
//#define SING_DEGENERATE_0121 121
//#define SING_DEGENERATE_0122 122
//#define SING_DEGENERATE_0221 221
//#define SING_DEGENERATE_0011 11
//#define SING_DEGENERATE_0110 110
//#define SING_DEGENERATE_0001 1
//#define SING_DEGENERATE_0111 111
#define SING_REGULAR 1
#define SING_CYCLIC 2
#define SING_SINGULAR 3
#define SING_DEGENERATE_0012 4
#define SING_DEGENERATE_0021 5
#define SING_DEGENERATE_0112 6
#define SING_DEGENERATE_0120 7
#define SING_DEGENERATE_0121 8
#define SING_DEGENERATE_0122 9
#define SING_DEGENERATE_0221 10
#define SING_DEGENERATE_0011 11
#define SING_DEGENERATE_0001 12
#define SING_DEGENERATE_0111 13
#define SING_CONSTANT 1
#define SING_IMPOSSIBLE 255

int gray_singularity_ordered(float a, float b, float c, float d)
{
	// a b
	// c d
	assert(a <= b); assert(a <= c); assert(a <= d);
	assert(b <= c);

	// 0 1      0 1      0 2
	// 2 3      3 2      3 1
	if (a < b && b < c && c < d) return SING_REGULAR;  // 0 1 2 3
	if (a < b && b < d && d < c) return SING_CYCLIC;   // 0 1 3 2
	if (a < d && d < b && b < c) return SING_SINGULAR; // 0 2 3 1

	// 0 0      0 0      0 1     0 1      0 1      0 1     0 2
	// 1 2      2 1      1 2     2 0      2 1      2 2     2 1
	if (a == b && b < c && c < d) return SING_DEGENERATE_0012;
	if (a == b && b < c && c > d) return SING_DEGENERATE_0021;
	if (a < b && b == c && c < d) return SING_DEGENERATE_0112;
	if (a == d && d < c && b < c) return SING_DEGENERATE_0120;
	if (a < b && b < c && d == b) return SING_DEGENERATE_0121;
	if (a < b && b < c && c == d) return SING_DEGENERATE_0122;
	if (a < b && b == c && c > d) return SING_DEGENERATE_0221;


	// 0 0      0 0      0 1
	// 1 1      0 1      1 1
	if (a == b && b < c && c == d) return SING_DEGENERATE_0011;
	if (a == b && b == c && c < d) return SING_DEGENERATE_0001;
	if (a < b && b == c && c == d) return SING_DEGENERATE_0111;
	if (a == b && b == c && c == d) return SING_CONSTANT; // 0 0 0 0
	fprintf(stderr, "impossible singularity: %g %g %g %g\n", a, b, c, d);
	return SING_IMPOSSIBLE;
}

int gray_singularity(float a, float b, float c, float d)
{
	if (a <= b && a <= c && a <= d) {
	       if (b <= c)
		       return gray_singularity_ordered(a, b, c, d);
	       else // swap diagonal
		       return gray_singularity(a, c, b, d);
	}
	else if (c <= a && c <= b && c <= d) // swap x axis
		return gray_singularity(c, d, a, b);
	else if (b <= a && b <= c && b <= d) // swap y axis
		return gray_singularity(b, a, d, c);
	else if (d <= a && d <= b && d <= c) // swap antidiagonal
		return gray_singularity(d, b, c, a);
	else
		assert(0);
}

void gray_singularities(char *y, float *x, int w, int h)
{
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		float a = x[w*(j+0)+(i+0)];
		float b = x[w*(j+0)+(i+1)];
		float c = x[w*(j+1)+(i+0)];
		float d = x[w*(j+1)+(i+1)];
		y[(w-1)*j+i] = gray_singularity(a, b, c, d);
	}
}

#include <stdlib.h>
#include <stdint.h>
#include "iio.h"

int main(int c, char **v)
{
	int w, h;
	float *x = iio_read_image_float("-", &w, &h);
	char *y = malloc(w*h);
	gray_singularities(y, x, w, h);
	iio_save_image_uint8_vec("-", (uint8_t*)y, w-1, h-1, 1);
	return 0;
}
