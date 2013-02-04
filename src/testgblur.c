// local lgblur (a section of the scale space)
// this version is optimized por generating multiple versions of the same image
// usage: lgblur2 infile.png varpat outpat FIRST LAST


#include <assert.h>
#include <math.h>
#include <stdio.h>


static float index_to_scale(int n, float s0, float sf, int i)
{
	float r = -1;
	if (i < 0) r = s0;
	else if (i >= n) r = sf;
	//else r = s0 * pow(sf/s0, i/(n-1.0));
	else r = s0 + i * (sf - s0)/(n - 1);
	assert(r >= s0);
	assert(r <= sf);
	return r;
}

static int scale_to_index(int n, float s0, float sf, float s)
{
	int r = -1;
	if (s <= s0) r = 0;
	else if (s >= sf) r = n - 1;
	//else r = lround((n - 1)*log(s/s0)/log(sf/s0));
	else r = lround((n - 1)*(s-s0)/(sf-s0));
	assert(r >= 0);
	assert(r < n);
	return r;
}

#define OMIT_GBLUR_MAIN
#include "gblur.c"

// utility headers used only in the "main" function
#include <stdio.h>

#include "iio.h"
#include "fail.c"
#include "xmalloc.c"

// fill a pyramid of images of the same size
void build_gaussian_pyramid(float *p, int npyr, float sfirst, float slast,
		float *x, int w, int h, int pd)
{
	assert(0 < sfirst);
	assert(sfirst < slast);
	for (int i = 0; i < npyr; i++)
	{
		float s = index_to_scale(npyr, sfirst, slast, i);
		fprintf(stderr, "s[%d] = %g\n", i, s);
		gblur(p+i*w*h*pd, x, w, h, pd, s);

		char buf[FILENAME_MAX];
		snprintf(buf, FILENAME_MAX, "/tmp/pyra_%03d.png", i);
		iio_save_image_float_vec(buf, p+i*w*h*pd, w, h, pd);
	}
}

static float imgdist_l1(float *a, float *b, int w, int h, int pd)
{
	long double r = 0;
	for (int i = 0; i < w*h*pd; i++)
		r += fabs(a[i] - b[i]);
	return r;
}

static float imgdist_l2(float *a, float *b, int w, int h, int pd)
{
	long double r = 0;
	for (int i = 0; i < w*h*pd; i++)
		r = hypot(r, a[i] - b[i]);
	return r;
}

static float imgdist_l3(float *a, float *b, int w, int h, int pd)
{
	long double r = 0;
	for (int i = 0; i < w*h*pd; i++)
		r += pow(fabs(a[i]-b[i]), 3);
	return pow(r, 1.0/3);
}


static void tryblur(float *sharp, float *blur, int w, int h, int pd, float s)
{
	float *bs = xmalloc(w*h*pd*sizeof*bs);
	gblur(bs, sharp, w, h, pd, s);
	float l1 = imgdist_l1(bs, blur, w, h, pd);
	float l2 = imgdist_l2(bs, blur, w, h, pd);
	float l3 = imgdist_l3(bs, blur, w, h, pd);
	fprintf(stderr, "%.9lf\tL1=%.9lf\tL2=%.9lf %.9lf\n", s, l1, l2, l3);
	free(bs);
}

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s sharp blur first last\n", *v);
		//                          0 1     2    3     4
		return 1;
	}
	char *filename_sharp = v[1];
	char *filename_blur = v[2];
	float first = atof(v[3]);
	float last = atof(v[4]);

	int w[2], h[2], pd[2];
	float *sharp = iio_read_image_float_vec(filename_sharp, w, h, pd);
	float *blur = iio_read_image_float_vec(filename_blur, w+1, h+1, pd+1);
	if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
		fail("input images size mismatch");

	int ntry = 10;
	for (int i = 0; i < ntry; i++)
	{
		float s = first + (i/(ntry-1.0)) * (last - first);
		tryblur(sharp, blur, *w, *h, *pd, s);
	}

	return 0;
}
