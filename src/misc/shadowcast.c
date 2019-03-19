#include <assert.h>
#include <math.h>
#include <stdio.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int insideP(int w, int h, int i, int j)
{
	return i>=0 && j>=0 && i<w && j<h;
}

static void canonical_bresenham_parkour(int *o, int w, int h, float p, float q)
{
	assert(q > 0);
	assert(q >= fabs(p));

	// fill bresenham offsets
	int t[h];
	for (int i = 0; i < h; i++)
		t[i] = lrint(i * p / q);

	// compute range of valid abscissae
	int i_min = 0;
	int i_max = w;
	if (p < 0) i_max += -t[h-1];
	if (p > 0) i_min -= t[h-1];

	// fill-in the paths
	int c = 0; // pixel counter
	for (int i = i_min; i < i_max; i++)
	{
		o[c++] = -1; // marks the beginning of a path
		for (int j = 0; j < h; j++)
			if (insideP(w, h, i+t[j], j))
				o[c++] = j*w + i+t[j];
	}

	assert(c < w*h + 2*(w+h));
}

static void compute_the_bresenham_parkour(
		int *o,    // output array of size w*h+2*(w+h), to be filled-in
		int w,     // width of the image domain
		int h,     // height of the image domain
		float p,   // cosine of direction
		float q    // sine of direction
		)
{
	int N = w*h + 2*(w+h);
	for (int i = 0; i < N; i++)
		o[i] = -1;

	if (fabs(p) <= q) // canonical case
		canonical_bresenham_parkour(o, w, h, p, q);
	else if (fabs(p) <= -q) { // swap the sign of q
		canonical_bresenham_parkour(o, w, h, p, -q);
		for (int i = 0; i < N; i++) if (o[i] >= 0)
			o[i] = (h - 1 - o[i]/w)*w + o[i]%w;
	} else if (fabs(q) < p) { // swap q and p
		canonical_bresenham_parkour(o, h, w, q, p);
		for (int i = 0; i < N; i++) if (o[i] >= 0)
			o[i] = (o[i]%h)*w + o[i]/h;
	} else if (fabs(q) < -p) { // swap q and -p
		//fprintf(stderr, "here pq = %g %g!", p, q);
		canonical_bresenham_parkour(o, h, w, q, -p);
		for (int i = 0; i < N; i++) if (o[i] >= 0)
			o[i] = (o[i]%h)*w + w - 1 - o[i]/h;
	}
}

#include "xmalloc.c"

static void cast_shadows(
		float *D,      // DEM raster data, to be filled-in with NAN
		int w,         // width of the raster
		int h,         // height of the raster
		float p,       // cosine of ground-projected sun direction
		float q,       // sine of ground-projected sun direction
		float a        // slope of the sun direction
		)
{
	fprintf(stderr, "cast shadows pqa = %g %g %g\n", p, q, a);
	int N = w*h + 2*(w+h);  // maximum possible length of parkour
	int *i = xmalloc(N*sizeof*i); // indices of Bresenham parkour
	compute_the_bresenham_parkour(i, w, h, p, q);
	int l = -1; // index of the first point on the current line
	for (int j = 0; j < N; j++)
	{
		if (i[j] < 0 && j+1<N && i[j+1] ) // mark the beginning of a line
			l = i[j+1];
		if (i[j] < 0 || l < 0) continue;
		assert(0 <= i[j]);
		assert(i[j] < w*h);
		if (l < 0) fprintf(stderr, "SHIT l=%d\n", l);
		assert(0 <= l);
		assert(l < w*h);
		float X = l % w;  // X of last occluding point
		float Y = l / w;  // Y of last occluding point
		float Z = D[l];   // Z of last occluding point
		float x = i[j] % w;  // x of current point
		float y = i[j] / w;  // y of current point
		float z = D[i[j]];   // z of current point
		if (Z - z > a * hypot(x-X, y-Y)) // occluded point
			D[i[j]] = NAN;
		else // new occluding point
			l = i[j];
	}
	free(i);
}

void cast_vertical_shadows(float *xx, int w, int h, float alpha)
{
	// pointer fo easy access ( x[j][i] == xx[j*w+i] )
	float (*x)[w] = (void*)xx;

	// compute slope of the rays
	float slope = tan(alpha * M_PI / 180);
	fprintf(stderr, "casting shadows with slope %g\n", slope);

	if (alpha <= 0)
		// process each column independently
		for (int i = 0; i < w; i++)
		{
			int l = 0;
			for (int j = 0; j < h; j++)
				if (x[j][i] < slope * (j - l) + x[l][i])
					x[j][i] = NAN;
				else
					l = j;
		}
	else
		for (int i = 0; i < w; i++)
		{
			float slop = tan((180-alpha)*M_PI/180);
			int l = h-1;
			for (int j = h-1; j >= 0; j--)
				if (x[j][i] < slop * (l - j) + x[l][i])
					x[j][i] = NAN;
				else
					l = j;
		}
}


//#define MAIN_VERTSHADOW
#define MAIN_THREEDSHADOW

#ifdef MAIN_VERTSHADOW
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"      // library for image input/output
#include "pickopt.c"  // function "pick_option" for processing args
int main(int c, char *v[])
{
	// process input arguments
	_Bool m = pick_option(&c, &v, "m", NULL);
	if (c < 2 || c > 4) {
		fprintf(stderr, "usage:\n\t%s alpha [dem_in [dem_out]]\n", *v);
		//                          0 1      2       3
		return 1;
	}
	float alpha = atof(v[1]);
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	if (pd != 1) {
		for (int i = 0; i < w*h; i++)
			for (int l = 1; l < pd; l++)
				x[i] += x[i+w*h*l];
		for (int i = 0; i < w*h; i++)
			x[i] /= pd;
	}

	// cast the vertical shadows
	cast_vertical_shadows(x, w, h, alpha);

	// if mask is requested, create a binary mask
	if (m) for (int i = 0; i < w*h; i++)
		x[i] = 255*isnan(x[i]);

	// save the output image
	iio_write_image_float_split(filename_out, x, w, h, 1);

	// cleanup (unnecessary) and exit
	return 0;
}
#endif//MAIN_VERTSHADOW


#ifdef MAIN_THREEDSHADOW
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"      // library for image input/output
int main(int c, char *v[])
{
	// process input arguments
	if (c < 4 || c > 6) {
		fprintf(stderr, "usage:\n\t%s p q r [dem_in [dem_out]]\n", *v);
		//                          0 1 2 3  4       5
		return 1;
	}
	float param_p = atof(v[1]);
	float param_q = atof(v[2]);
	float param_r = atof(v[3]);
	char *filename_in  = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	// read input dem
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	if (pd != 1) {
		for (int i = 0; i < w*h; i++)
			for (int l = 1; l < pd; l++)
				x[i] += x[i+w*h*l];
		for (int i = 0; i < w*h; i++)
			x[i] /= pd;
	}


	cast_shadows(x, w, h, param_p, param_q, param_r);


	//int *o = xmalloc( (w*h + 2*(w+h)) * sizeof*o);
	//compute_the_bresenham_parkour(o, w, h, param_p, param_q);
	//for (int i = 0; i < w*h; i++)
	//	x[i] = -1;
	//int cx = 0;
	//for (int i = 0; i < w*h+2*(w+h); i++)
	//	if (o[i] >= 0)
	//	{
	//		assert(o[i] < w*h);
	//		x[o[i]] = cx++;
	//	}
	//free(o);

	//// cast the vertical shadows
	//cast_vertical_shadows(x, w, h, alpha);

	//// if mask is requested, create a binary mask
	//if (m) for (int i = 0; i < w*h; i++)
	//	x[i] = 255*isnan(x[i]);

	// save the output image
	iio_write_image_float_split(filename_out, x, w, h, 1);

	// cleanup (unnecessary) and exit
	free(x);
	return 0;
}
#endif//MAIN_THREEDSHADOW


