// this program dumps the given image into the terminal, using
// utf-8 characters colored with 24 bit capabilities
//
// no resizing of the input image is done, your image must have
// a width smaller than 80
//
// EXAMPLE:
// zoombil 80 -1 /tmp/lenac.png - | qauto -i | ./bin/idump


#include <stdint.h>
#include <stdio.h>
#include "iio.h"

static void idump_pair(uint8_t *a, uint8_t *b)
{
	printf("\x1b[38;2;%d;%d;%dm\x1b[48;2;%d;%d;%dm▀",
			a[0], a[1], a[2], b[0], b[1], b[2]);

// https://en.wikipedia.org/wiki/ANSI_escape_code#24-bit
//
//   ESC[ … 38;2;<r>;<g>;<b> … m Select RGB foreground color
//   ESC[ … 48;2;<r>;<g>;<b> … m Select RGB background color
//
// NOTE: \x1b == ESC
}


static void idump(uint8_t *x, int w, int h)
{
	for (int j = 0; j < h; j += 2)
	{
		for (int i = 0; i < w; i++)
		{
			uint8_t *a = x + 3 * (w*(j+0) + i);
			uint8_t *b = x + 3 * (w*(j+1) + i);
			idump_pair(a, b);
		}
		printf("\x1b[0m\n");
	}
}

int main(int c, char *v[])
{
	int w, h, pd;
	uint8_t *x = iio_read_image_uint8_vec("-", &w, &h, &pd);
	if (pd != 3) return 1;
	idump(x, w, h);
	return 0;
}
