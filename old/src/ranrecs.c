#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"



#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)
#define FORK(n) for(int k=0;k<(n);k++)


#define BAD_MIN(a,b) (b)<(a)?(b):(a)
#define BAD_MAX(a,b) (b)>(a)?(b):(a)

#include "fragments.c"

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

int bound(int a, int x, int b)
{
	if (b < a) return bound(b, a, x);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

int main(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t%s w h nf nrecs outpat seed\n", *v);
		//                          0 1 2 3  4     5      6
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int nf = atoi(v[3]);
	int nrecs = atoi(v[4]);
	char *outpat = v[5];
	int seed = atoi(v[6]);
	srand(seed);

	int rec[nrecs][9]; // (x0, y0, xf, yf, r, g, b, vx, vy)
	                   //  0   1   2   3   4  5  6  7   8
	FORL(nrecs) {
		rec[l][0] = randombounds(-10, w+10);
		rec[l][1] = randombounds(-10, h+10);
		rec[l][2] = rec[l][0] + randombounds(20, w/2);
		rec[l][3] = rec[l][1] + randombounds(20, h/2);
		rec[l][4] = randombounds(0, 255);
		rec[l][5] = randombounds(0, 255);
		rec[l][6] = randombounds(0, 255);
		rec[l][7] = randombounds(-5, 5);
		rec[l][8] = randombounds(-5, 5);
	}

	char (*x)[w][3] = xmalloc(w*h*3);
	FORK(nf) {
		FORJ(h) FORI(w) FORL(3) x[j][i][l] = 0;
		FORL(nrecs) {
			int r[4];
			r[0] = rec[l][0] + k*rec[l][7];
			r[1] = rec[l][1] + k*rec[l][8];
			r[2] = rec[l][2] + k*rec[l][7];
			r[3] = rec[l][3] + k*rec[l][8];
			r[0] = bound(0, r[0], w-1);
			r[1] = bound(0, r[1], h-1);
			r[2] = bound(0, r[2], w-1);
			r[3] = bound(0, r[3], h-1);
			for (int i = r[0]; i <= r[2]; i++)
			for (int j = r[1]; j <= r[3]; j++)
			for (int s = 0; s < 3; s++)
				x[j][i][s] = rec[l][4+s];
		}
		char buf[0x100];
		snprintf(buf, 0x100, outpat, k);
		iio_write_image_uint8_vec(buf, (void*)x, w, h, 3);
	}
	xfree(x);

	return 0;
}
