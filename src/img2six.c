#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "iio.h"

//// generic 8-bit RGB palette (3+3+2 bits)
//static void fill_palette(uint8_t p[0x100][3])
//{
//	for (int i = 0; i < 0x100; i++)
//	{
//		int r = i/32;
//		int g = (i/4) % 8;
//		int b = i % 4;
//		if( i != 32*r + 4*g + b)
//			printf("i=%d rgb = %d %d %d\n", i, r, g, b);
//		assert( i == 32*r + 4*g + b);
//		//int R = floor(25.0 * r / 64.0);
//		//int G = floor(25.0 * g / 64.0);
//		//int B = floor(25.0 * b / 64.0);
//		//assert(R >= 0); assert(R < 100);
//		//assert(G >= 0); assert(G < 100);
//		//assert(B >= 0); assert(B < 100);
//		p[i][0] = r;// * 8;
//		p[i][1] = g;// * 8;
//		p[i][2] = b;// * 4;
//	}
//}

static int pal_idx(uint8_t *rgb)
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

	//uint8_t p[0x100][3];
	//fill_palette(p);

	// beginning of sixel block
	printf("\ePq\n");
	//printf("<ESC>Pq\n");

	// dump palette (RGB 3+3+2 bits palette)
	for (int i = 0; i < 0x100; i++)
	{
		int r = i/32;
		int g = (i/4) % 8;
		int b = i % 4;
		printf("#%d;2;%d;%d;%d", i, //i/32, (i/4)%8, i%4);
				(int)floor(100.0*r/7.0),
				(int)floor(100.0*g/7.0),
				(int)floor(100.0*b/3.0)
		);
	}

	// dump sixels
	for (int jj = 0; jj < h/6; jj++) // traverse sixel stripes
	{
		// find colors that appear on this stripe
		int m[0x100] = {0}, cnt = 0;
		for (int i = 0; i < 6*w; i++)
		{
			int k = pal_idx(x+3*(6*jj*w+i));
			assert(k >= 0); assert(k < 256);
			if (!m[k]) cnt += 1;
			m[k] += 1;
		}

		// for each color that appears, dump its binary mask
		for (int k = 0; k < 0x100; k++)
		if (m[k])
		{
			uint8_t buf[w]; // sixels before run-length compression
			cnt -= 1;
			for (int i = 0; i < w; i++) // traverse sixels
			{
				int s[6]; // sixel binary mask
				for (int l = 0; l < 6; l++)
					s[l] = k==pal_idx(x+3*((6*jj+l)*w+i));
				int S = s[0] + 2*( // sixel index
					s[1] + 2*(
					s[2] + 2*(
					s[3] + 2*(
					s[4] + 2*( s[5] )))));
				assert(S >= 0 && S < 64);
				buf[i] = S + 63;
				//fprintf(stderr, "buf[i=%d]=%d\n", i, buf[i]);
				//printf("%c", S + 63); // dump sixel
			}

			// run-lenght encode
			uint8_t rle[w];
			int rep[w] , idx = 0;
			for (int i = 0; i < w; i++)
				rep[i] = 1;
			rle[0] = *buf;
			for (int i = 1; i < w; i++)
				if (buf[i] == rle[idx])
					rep[idx] += 1;
				else
					rle[++idx] = buf[i];

			//fprintf(stderr, "w=%d idx=%d\n", w, idx);
			//fprintf(stderr, "buf = %d %d %d ...\n", buf[0], buf[1], buf[2]);
			//fprintf(stderr, "rle = %d %d %d ...\n", rle[0], rle[1], rle[2]);
			//fprintf(stderr, "rep = %d %d %d ...\n", rep[0], rep[1], rep[2]);

			// dump sixels
			printf("#%d", k); // select color k
			//for (int i = 0; i < w; i++)
			//	printf("%c", buf[i]);
			for (int i = 0; i <= idx; i++)
				if (rep[idx] < 3)
					for (int k = 0; k < rep[i]; k++)
						printf("%c", rle[i]);
				else
					printf("!%d%c", rep[i], rle[i]);
			printf(cnt ? "$\n" : "-\n"); // signal end of line
		}
	}

	// end of sixel block
	printf("\e\\");
	//printf("<ESC>\\");
}

static int pidx(uint8_t *rgb)
{
	return 32*(rgb[0]/32) + 4*(rgb[1]/32) + rgb[2]/64;
}

void dump_sixels_to_stdout_rgb2(uint8_t *x, int w, int h)
{
	// beginning of sixel block
	printf("\ePq\n");

	// dump palette (RGB 3+3+2 bits palette)
	for (int i = 0; i < 0x100; i++)
		printf("#%d;2;%d;%d;%d", i,
				(int)(14.2857*(i/32)),
				(int)(14.2857*((i/4)%8)),
				(int)(33.3333*(i%4))
		);

	// dump sixels
	for (int jj = 0; jj < h/6; jj++) // traverse sixel stripes
	{
		// find colors that appear on this stripe
		int m[0x100] = {0}, cnt = 0;
		for (int i = 0; i < 6*w; i++)
		{
			int k = pidx(x+3*(6*jj*w+i));
			if (!m[k]) cnt += 1;
			m[k] += 1;
		}

		// for each color that appears, dump its binary mask
		for (int k = 0; k < 0x100; k++)
		if (m[k])
		{
			uint8_t buf[w]; // sixels before run-length compression
			cnt -= 1;
			for (int i = 0; i < w; i++) // traverse sixels
			{
				int s[6]; // sixel binary mask
				for (int l = 0; l < 6; l++)
					s[l] = k == pidx(x+3*((6*jj+l)*w+i));
				int S = s[0] + 2*( // sixel index
					s[1] + 2*(
					s[2] + 2*(
					s[3] + 2*(
					s[4] + 2*( s[5] )))));
				buf[i] = S + 63;
			}

			// run-lenght encode
			uint8_t rle[w];
			int rep[w] , idx = 0;
			for (int i = 0; i < w; i++)
				rep[i] = 1;
			rle[0] = *buf;
			for (int i = 1; i < w; i++)
				if (buf[i] == rle[idx])
					rep[idx] += 1;
				else
					rle[++idx] = buf[i];

			// dump sixels
			printf("#%d", k); // select color k
			for (int i = 0; i <= idx; i++)
				if (rep[idx] < 3)
					for (int k = 0; k < rep[i]; k++)
						printf("%c", rle[i]);
				else
					printf("!%d%c", rep[i], rle[i]);
			printf(cnt ? "$\n" : "-\n"); // signal end of line
		}
	}

	// end of sixel block
	printf("\e\\");
}

void dump_sixels_to_stdout_rgb3(uint8_t *x, int w, int h)
{
	printf("\ePq\n");
	for (int i = 0; i < 0x100; i++)
		printf("#%d;2;%d;%d;%d", i,
				(int)(14.2857*(i/32)),
				(int)(14.2857*((i/4)%8)),
				(int)(33.3333*(i%4))
		);
	for (int jj = 0; jj < h/6; jj++)
	{
		int m[0x100] = {0}, cnt = 0;
		for (int i = 0; i < 6*w; i++)
		{
			int k = pidx(x+3*(6*jj*w+i));
			if (!m[k]) cnt += 1;
			m[k] += 1;
		}
		for (int k = 0; k < 0x100; k++)
		if (m[k])
		{
			uint8_t buf[w];
			cnt -= 1;
			for (int i = 0; i < w; i++)
			{
				int s[6];
				for (int l = 0; l < 6; l++)
					s[l] = k == pidx(x+3*((6*jj+l)*w+i));
				int S = s[0] + 2*( s[1] + 2*( s[2] + 2*(
					s[3] + 2*( s[4] + 2*( s[5] )))));
				buf[i] = S + 63;
			}
			uint8_t rle[w];
			int rep[w] , idx = 0;
			for (int i = 0; i < w; i++)
				rep[i] = 1;
			rle[0] = *buf;
			for (int i = 1; i < w; i++)
				if (buf[i] == rle[idx])
					rep[idx] += 1;
				else
					rle[++idx] = buf[i];
			printf("#%d", k);
			for (int i = 0; i <= idx; i++)
				if (rep[idx] < 3)
					for (int k = 0; k < rep[i]; k++)
						printf("%c", rle[i]);
				else
					printf("!%d%c", rep[i], rle[i]);
			printf(cnt ? "$\n" : "-\n");
		}
	}
	printf("\e\\");
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
		dump_sixels_to_stdout_rgb2(x, w, h);
	else
		return fprintf(stderr, "d = %d not supported yet\n", d);

	return 0;
}
