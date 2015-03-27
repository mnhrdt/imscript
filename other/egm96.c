#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

// location of data file
#define EGM96_025_TIF "/home/coco/.srtm4/egm96_025.tiff"


// data loaded in RAM
static float global_egm96_025_data[721*1440] = { 0 };

// function to evaluate the bilinear interpolation inside a 2x2 cell
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == (unsigned int)p)
			r = 0;
	}
	return r;
}

// getpixel with periodic boundary condtitions
static float getpixel_per(float *x, int w, int h, int i, int j)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	assert(i >= 0 && i < w);
	assert(j >= 0 && j < h);
	return x[j*w+i];
}

// bilinear interpolation at the given sub-pixel position
static float bilinear_interpolation_at(float *x, int w, int h, float p, float q)
{
	if (!x) return NAN;
	int ip = floor(p);
	int iq = floor(q);
	float a = getpixel_per(x, w, h, ip  , iq  );
	float b = getpixel_per(x, w, h, ip+1, iq  );
	float c = getpixel_per(x, w, h, ip  , iq+1);
	float d = getpixel_per(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// evaluate the egm96 geoid at the requested site, expressed in degrees
double egm96(double longitude, double latitude)
{
	static bool firstcall = true;
	int w = 1440;
	int h = 720;
	if (firstcall) {
		int ww, hh;
		float *tmp = iio_read_image_float(EGM96_025_TIF, &ww, &hh);
		if (!tmp || ww != 1+w || hh != 1+h)
			fprintf(stderr, "WARNING: could not read EGM96 data\n");
		else if (tmp) {
			// remove repeated last column
			float *x = global_egm96_025_data;
			for (int j = 0; j < h; j++)
			for (int i = 0; i < w; i++)
				x[j*w + i] = tmp[j*(w + 1) + i];
			free(tmp);
		}
		firstcall = false;
	}
	float fi = 4 * longitude;
	float fj = 4 * (90 - latitude);
	return bilinear_interpolation_at(global_egm96_025_data, w, h, fi, fj);
}

#ifdef MAIN_EGM96
int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s longitude latitude\n", *v);
		return 1;
	}
	double lon = atof(v[1]);
	double lat = atof(v[2]);
	double r = egm96(lon, lat);
	printf("%g\n", r);
	return 0;
}
#endif//MAIN_EGM96
