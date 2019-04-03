#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "iio.h"

// generic 8-bit RGB palette (3+3+2 bits)
static void fill_palette(uint8_t p[0x100][3])
{
	for (int i = 0; i < 0x100; i++)
	{
		int r = i/32;
		int g = (i/4) % 8;
		int b = i % 4;
		if( i != 32*r + 4*g + b)
			printf("i=%d rgb = %d %d %d\n", i, r, g, b);
		assert( i == 32*r + 4*g + b);
		//int R = floor(25.0 * r / 64.0);
		//int G = floor(25.0 * g / 64.0);
		//int B = floor(25.0 * b / 64.0);
		//assert(R >= 0); assert(R < 100);
		//assert(G >= 0); assert(G < 100);
		//assert(B >= 0); assert(B < 100);
		p[i][0] = r;// * 8;
		p[i][1] = g;// * 8;
		p[i][2] = b;// * 4;
	}
}

static int get_pal_idx(uint8_t *rgb)
{
	uint8_t r = rgb[0]/32;
	uint8_t g = rgb[1]/32;
	uint8_t b = rgb[2]/64;
	int i = 32*r + 4*g + b;
	if (i > 255)
		printf("BAD i=%d rgb=%d %d %d\n", i, r, g, b);
	return i;
}

void dump_sixels_to_stdout_rgb(uint8_t *x, int w, int h)
{
	//uint8_t *debug = malloc(3*w*h);
	//for (int i = 0; i < 3*w*h; i++) debug[i] = 0;

	uint8_t p[0x100][3];
	fill_palette(p);

	// beginning of sixel block
	printf("\ePq\n");
	//printf("<ESC>Pq\n");

	// dump palette
	for (int i = 0; i < 0x100; i++)
		printf("#%d;2;%d;%d;%d", i,
				(int)floor(100.0*p[i][0]/7.0),
				(int)floor(100.0*p[i][1]/7.0),
				(int)floor(100.0*p[i][2]/3.0)
		      );

	// dump sixels
	for (int jj = 0; jj < h/6; jj++) // traverse sixel stripes
	{
		// find colors that appear on this stripe
		int m[0x100] = {0}, cnt = 0;
		for (int i = 0; i < 6*w; i++)
		{
			int k = get_pal_idx(x+3*(6*jj*w+i));
			assert(k >= 0); assert(k < 256);
			if (!m[k]) cnt += 1;
			m[k] += 1;
		}

		// for each color that appears, dump its binary mask
		for (int k = 0; k < 0x100; k++)
		if (m[k])
		{
			cnt -= 1;
			printf("#%d", k); // select color k
			for (int i = 0; i < w; i++) // traverse sixels
			{
				int s[6]; // sixel binary mask
				for (int l = 0; l < 6; l++)
					s[l] = k==get_pal_idx(x+3*((6*jj+l)*w+i));
				int S = s[0] + 2*( // sixel index
					s[1] + 2*(
					s[2] + 2*(
					s[3] + 2*(
					s[4] + 2*( s[5] )))));
				assert(S >= 0 && S < 64);
				printf("%c", S + 63); // dump sixel
			}
			printf(cnt ? "$\n" : "-\n"); // signal end of line
		}
	}

	// end of sixel block
	printf("\e\\");
	//printf("<ESC>\\");
}

int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s img.png\n", *v);
	//                                         0 1
	char *filename_in = v[1];

	int w, h, d;
	uint8_t *x = iio_read_image_uint8_vec(filename_in, &w, &h, &d);

	if (d == 3)
		dump_sixels_to_stdout_rgb(x, w, h);
	else
		return fprintf(stderr, "d = %d not supported yet\n", d);

	return 0;
}
