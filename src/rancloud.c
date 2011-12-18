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

static void raninuball(float *x)
{
	float r;
	do {
		FORI(3) x[i] = 2*random_uniform() - 1;
		r = hypot(x[0], hypot(x[1],x[2]));
	} while(r > 1);
}

int main(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t%s w h nf ndots outpat seed\n", *v);
		//                          0 1 2 3  4     5      6
		return 1;
	}
	int w = atoi(v[1]);
	//int h = atoi(v[2]);
	int nf = atoi(v[3]);
	int ndots = atoi(v[4]);
	char *outpat = v[5];
	int seed = atoi(v[6]);
	srand(seed);

	float dot[ndots][3];
	FORL(ndots) raninuball(dot[l]);

	unsigned char (*x)[w] = xmalloc(w*w);
	FORK(nf) {
		float theta = k*M_PI/(0.5*nf);
		float vec[2] = {cos(theta), sin(theta)};
		FORJ(w) FORI(w) x[j][i] = 255;
		FORL(ndots) {
			float X = dot[l][0]*vec[0] + dot[l][1]*vec[1];
			float Y = dot[l][2];
			int iX = (w - 1) * (X + 1) / 2;
			int iY = (w - 1) * (Y + 1) / 2;
			iX = bound(0, iX, w-1);
			iY = bound(0, iY, w-1);
			x[iY][iX] = 0;
		}
		char buf[0x100];
		snprintf(buf, 0x100, outpat, k);
		iio_save_image_uint8_vec(buf, (void*)x, w, w, 1);
	}
	xfree(x);

	return 0;
}
