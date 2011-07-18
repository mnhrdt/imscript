#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"

static float getsample(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

static void setsample(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

static void setpixel(float *x, int w, int h, int i, int j, float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return;
	x[i + j*w] = v;
}

static float evaluate_bilinear(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static void bilinear_zoom(float *X, int W, int H, float *x, int w, int h)
{
	float wfactor = w/(float)W;
	float hfactor = h/(float)H;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		float p = i*wfactor;
		float q = j*hfactor;
		int ip = p;
		int iq = q;
		float a = getpixel(x, w, h, ip  , iq  );
		float b = getpixel(x, w, h, ip+1, iq  );
		float c = getpixel(x, w, h, ip  , iq+1);
		float d = getpixel(x, w, h, ip+1, iq+1);
		float r = evaluate_bilinear(a, b, c, d, p-ip, q-iq);
		setpixel(X, W, H, i, j, r);
	}
}


//static float *bilinear_zoom(float *x, int w, int h, int pd, int ow, int oh)
//{
//	float *y = malloc(ow*oh*pd*sizeof*y);
//	float wfactor = w/(float)ow;
//	float hfactor = h/(float)oh;
//	for (int j = 0; j < oh; j++)
//	for (int i = 0; i < ow; i++)
//	for (int l = 0; l < pd; l++)
//	{
//		float p = i*wfactor;
//		float q = j*hfactor;
//		int ip = p;
//		int iq = q;
//		float a = getsample(x, w, h, pd, ip  , iq  , l);
//		float b = getsample(x, w, h, pd, ip+1, iq  , l);
//		float c = getsample(x, w, h, pd, ip  , iq+1, l);
//		float d = getsample(x, w, h, pd, ip+1, iq+1, l);
//		float r = evaluate_bilinear(a, b, c, d, p-ip, q-iq);
//		setsample(y, ow, oh, pd, i, j, l, r);
//	}
//	return y;
//}

static void bilinear_interpolation(int pd, float (**y)[pd], int ow, int oh,
		float (**x)[pd], int w, int h)
{
	float wfactor = w/(float)ow;
	float hfactor = h/(float)oh;
	float a, b, c, d;
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < pd; l++)
	{
		float p = i*wfactor;
		float q = j*hfactor;
		int ip = p;
		int iq = q;
		//fprintf(stderr, "i j p q ip iq %d %d %g %g %d %d\n", i, j, p, q, ip, iq);
		a = x[iq][ip][l];
		if (ip+1 < w) b = x[iq][ip+1][l];
		if (iq+1 < h) c = x[iq+1][ip][l];
		if (ip+1 < w && iq+1 < h) d = x[iq+1][ip+1][l];
		float r = evaluate_bilinear(a, b, c, d, p-ip, q-iq);
		y[j][i][l] = r;
	}
}

//// input: n vectors of dimension d and component size s
//// output (in place): d arrays of length n and element size s
//static void act_of_violence(void *data, int n, int d, int s)
//{
//	void *tmp = malloc(n*d*s);
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < d; j++)
//		{
//			void *from = s*(d*i + j) + (char *)data;
//			void *to =   s*(n*j + i) + (char *)tmp;
//			memcpy(to, from, s);
//		}
//
//	memcpy(data, tmp, n*d*s);
//	free(tmp);
//}

//int main_flat(int c, char *v[])
//{
//	if (c != 3) {
//		fprintf(stderr, "usage:\n\t%s width height < in > out\n", *v);
//		return EXIT_FAILURE;
//	}
//	int ow = atoi(v[1]);
//	int oh = atoi(v[2]);
//
//	int w, h, pd;
//	float *x = iio_read_image_float_vec("-", &w, &h, &pd);
//	float *y = bilinear_zoom(x, w, h, pd, ow, oh);
//
//	iio_save_image_float_vec("-", y, ow, oh, pd);
//
//	free(y);
//	free(x);
//
//	return EXIT_SUCCESS;
//}

// alloc a 2d matrix contiguously (wxh elements of size n)
static void *matrix_build(int w, int h, size_t n)
{
	size_t p = sizeof(void *);
	char *r = malloc(h*p + w*h*n);
	for (int i = 0; i < h; i++)
		*(void **)(r + i*p) = r + h*p + i*w*n;
	return r;
}

static int main_gray(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s width height in out\n", *v);
		return EXIT_FAILURE;
	}

	int W = atoi(v[1]);
	int H = atoi(v[2]);

	int w, h;
	float *x = iio_read_image_float(v[3], &w, &h);
	iio_save_image_float("/tmp/caca.png", x, w, h);
	float *y = malloc(W*H*sizeof*y);
	bilinear_zoom(y, W, H, x, w, h);

	iio_save_image_float(v[4], y, W, H);

	free(y);
	free(x);

	return EXIT_SUCCESS;
}

static int main_vec(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s width height in out\n", *v);
		return EXIT_FAILURE;
	}
	int ow = atoi(v[1]);
	int oh = atoi(v[2]);

	int w, h, pd;
	void *data = iio_read_image_float_matrix_vec(v[3], &w, &h, &pd);
	float (**x)[pd] = data;


	float (**y)[pd] = matrix_build(ow, oh, sizeof**y);
	bilinear_interpolation(pd, y, ow, oh, x, w, h);

	iio_save_image_float_vec(v[4], y[0][0], ow, oh, pd);

	free(x);
	free(y);

	return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
	return main_vec(c, v);
	/*
	return main_gray(c, v);
	*/
}
