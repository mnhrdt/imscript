#include "xmalloc.c"

typedef float (*getpixel_operator)(float*,int,int,int,int);

static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

static float dx(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return p(x,w,h, i+1, j) - p(x,w,h, i, j);
}

static float dy(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return p(x,w,h, i, j+1) - p(x,w,h, i, j);
}

static float cdx(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return p(x,w,h, i+1, j) - p(x,w,h, i-1, j);
}

static float cdy(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return p(x,w,h, i, j+1) - p(x,w,h, i, j-1);
}

static float sobel_dx(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return 0
	+ 2*p(x,w,h, i+1, j) - 2*p(x,w,h, i-1, j)
	+ p(x,w,h, i+1, j+1) - p(x,w,h, i-1, j+1)
	+ p(x,w,h, i+1, j-1) - p(x,w,h, i-1, j-1);
}

static float sobel_dy(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;
	return 0
	+ 2*p(x,w,h, i, j+1) - 2*p(x,w,h, i, j-1)
	+ p(x,w,h, i+1, j+1) - p(x,w,h, i+1, j-1)
	+ p(x,w,h, i-1, j+1) - p(x,w,h, i-1, j-1);
}

static void vector_product(float ab[3], float a[3], float b[3])
{
	// i  j  k
	// a0 a1 a2
	// b0 b1 b2
	ab[0] = a[1]*b[2] - b[1]*a[2];
	ab[1] = a[2]*b[0] - b[2]*a[0];
	ab[2] = a[0]*b[1] - b[0]*a[1];
}

static float scalar_product(float *a, float *b, int n)
{
	return n ? scalar_product(a+1, b+1, n-1) + *a * *b: 0;
}

static void get_normal_at(float n[3], float *x, int w, int h, int i, int  j)
{
	float vdx[3] = {1, 0, cdx(x,w,h,i,j)};
	float vdy[3] = {0, 1, cdy(x,w,h,i,j)};
	vector_product(n, vdx, vdy);
}

void hshading(float *y, float *x, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < h; i++)
	{
		float normal[3];
		float sun[3] = {-1, -1, 1};
		get_normal_at(normal, x, w, h, i, j);
		float sp = scalar_product(normal, sun, 3);
		y[j*w + i] = sp;
	}
}


#include "iio.h"
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s height.tiff view.png\n", *v);
		//                         0  1           2
		return 1;
	}
	char *filename_in = v[1];
	char *filename_out = v[2];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	fprintf(stderr, "x w h = %p %d %d\n", (void*)x, w, h);

	float *y = xmalloc(w*h*sizeof*y);

	hshading(y, x, w, h);

	iio_write_image_float_vec(filename_out, y, w, h, 1);

	free(x);
	free(y);
	return 0;
} 
