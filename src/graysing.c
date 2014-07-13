#include <assert.h>
#include <stdio.h>


// 0 1
// 2 3
//
// 0 1
// 3 2
//
// 0 2
// 3 1
//
// 0 0
// 1 2
//
// 0 1
// 1 2
//
// 0 1
// 2 2

#define SING_REGULAR 0
#define SING_CYCLIC 1
#define SING_SINGULAR 2
#define SING_DEGENERATE 3
#define SING_DEGENERATE_SINGULAR 4
#define SING_DOUBLY_DEGENERATE 5
#define SING_DOUBLY_DEGENERATE_SINGULAR 6
#define SING_TRIPLY_DEGENERATE_LOW 7
#define SING_TRIPLY_DEGENERATE_HIGH 8
#define SING_CONSTANT 9
#define SING_IMPOSSIBLE 10

// a b
// c d


int gray_singularity_ordered(float a, float b, float c, float d)
{
	assert(a <= b); assert(a <= c); assert(a <= d);
	assert(b <= c);

	if (a < b && b < c && c < d) return SING_REGULAR;
	if (a < b && b < d && d < c) return SING_CYCLIC;
	if (a < d && d < b && b < c) return SING_SINGULAR;
	if (a == b && b < c && c < d) return SING_DEGENERATE;
	if (a == d && d < c && c < b) return SING_DEGENERATE_SINGULAR;
	if (a == b && b < c && c == d) return SING_DOUBLY_DEGENERATE;
	if (a == d && d < c && c == b) return SING_DOUBLY_DEGENERATE_SINGULAR;
	if (a == b && b == c && c < d) return SING_TRIPLY_DEGENERATE_LOW;
	if (a < b && b == c && c == d) return SING_TRIPLY_DEGENERATE_HIGH;
	if (a == b && b == c && c == d) return SING_CONSTANT;
	fprintf(stderr, "%g %g %g %g\n", a, b, c, d);
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
	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		float a = x[w*(j+0)+(i+0)];
		float b = x[w*(j+0)+(i+1)];
		float c = x[w*(j+1)+(i+0)];
		float d = x[w*(j+1)+(i+1)];
		y[w*j+i] = gray_singularity(a, b, c, d);
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
	iio_save_image_uint8_vec("-", (uint8_t*)y, w, h, 1);
	return 0;
}
