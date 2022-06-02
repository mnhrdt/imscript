#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

// tmp dbg
#include "iio.h"      // library for image input/output


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//static int insideP(int w, int h, int i, int j)
//{
//	return i>=0 && j>=0 && i<w && j<h;
//}
#define INSIDEP(w,h,i,j) (i>=0 && j>=0 && i<w && j<h)

//static void fill_canonical_bresenham_seeds(
//		float opq[2]              // (output) traversal step
//		int *oxy,                 // (output) coordinates of the seeds
//		int w, int h,             // width, height of the image domain
//		float p, float q          // normalized direction
//		)
//{
//	assert(q > 0);
//	assert(q >= fabs(p));
//
//	opq[0] = p/q;  // x-step (common for all paths)
//	opq[1] = 1;    // y-step (common for all paths)
//
//	// fill bresenham offsets
//	int t[h];
//	for (int i = 0; i < h; i++)
//		t[i] = lrint(i * p / q);
//
//	// compute range of valid abscissae
//	int i_min = 0;
//	int i_max = w;
//	if (p < 0) i_max += -t[h-1];
//	if (p > 0) i_min -= t[h-1];
//
//	// fill-in the beginning of each path
//	int c = 0; // pixel counter
//	for (int i = i_min; i < i_max; i++)
//	{
//		for (int j = 0; j < h; j++)
//			if (INSIDEP(w, h, i+t[j], j))
//			{
//				o[c++] = j*w + i+t[j];
//			}
//	}
//
//	assert(c < w + h);
//}
// sorry, I can't do this right now

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
			if (INSIDEP(w, h, i+t[j], j))
				o[c++] = j*w + i+t[j];
	}

	assert(c < w*h + 2*(w+h));
}

static void canonical_bresenham_parkour_hack(
		int *o,            // output pixel indices
		float *oxy,        // output actual coordinates
		int w, int h,      // image domain
		float p, float q   // direction
		)
{
	assert(q > 0);
	assert(q >= fabs(p));

	// fill bresenham offsets
	int t[h];
	for (int j = 0; j < h; j++)
		t[j] = lrint(j * p / q);
		//t[j] = floor(j * p / q);

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
			if (INSIDEP(w, h, i+t[j], j))
			{
				oxy[2*c+0] = i + j*p/q;
				oxy[2*c+1] = j;
				o[c] = j*w + i+t[j];
				c += 1;
			}
	}

	assert(c < w*h + 2*(w+h));
}

static void compute_the_bresenham_parkour_hack(
		int *o,     // output array of size w*h+2*(w+h), to be filled-in
		float *oxy, // output actual coordinates
		int w,      // width of the image domain
		int h,      // height of the image domain
		float p,    // cosine of direction
		float q     // sine of direction
		)
{
	int N = w*h + 2*(w+h);
	for (int i = 0; i < N; i++)
		o[i] = -1;

	if (fabs(p) <= q) // canonical case
	{
		canonical_bresenham_parkour_hack(o, oxy, w, h, p, q);
		fprintf(stderr, "SWAP nothing\n");
	}
	else if (fabs(p) <= -q) { // swap the sign of q
		canonical_bresenham_parkour_hack(o, oxy, w, h, p, -q);
		fprintf(stderr, "SWAP the sign of q\n");
		for (int i = 0; i < N; i++)
			if (o[i] >= 0) {
				o[i] = (h - 1 - o[i]/w)*w + o[i]%w;
				oxy[2*i+1] = h - 1 - oxy[2*i+1];
			}
	} else if (fabs(q) < p) { // swap q and p
		canonical_bresenham_parkour_hack(o, oxy, h, w, q, p);
		fprintf(stderr, "SWAP q and p\n");
		for (int i = 0; i < N; i++)
			if (o[i] >= 0) {
				o[i] = (o[i]%h)*w + o[i]/h;
				float t = oxy[2*i+0];
				oxy[2*i+0] = oxy[2*i+1];
				oxy[2*i+1] = t;
			}
	} else if (fabs(q) < -p) { // swap q and -p
		//fprintf(stderr, "here pq = %g %g!", p, q);
		canonical_bresenham_parkour_hack(o, oxy, h, w, q, -p);
		fprintf(stderr, "SWAP q and -p\n");
		for (int i = 0; i < N; i++)
			if (o[i] >= 0) {
				o[i] = (o[i]%h)*w + w - 1 - o[i]/h;
				float t = oxy[2*i+1];
				oxy[2*i+1] = oxy[2*i+0];
				oxy[2*i+0] = w - 1 - t;
			}
	}

	// sanity check
	for (int i = 0; i < N; i++)
		if (o[i] >= 0)
		{
			int ix = o[i] % w;
			int iy = o[i] / w;
			assert(0 <= ix && ix < w);
			assert(0 <= iy && iy < h);
			float x = oxy[2*i+0];
			float y = oxy[2*i+1];
			if (fabs(x - ix) > 0.5 || fabs(y-iy) > 0.5) fprintf(stderr, "i=%d, o[i]=%d, ixy=%d %d xy=%g %g\n", i, o[i], ix, iy, x, y);
			assert(fabs(x - ix) < 0.6);
			assert(fabs(y - iy) < 0.6);
		}
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

//static void compute_the_bresenham_seeds(
//		int *o,    // output array of size w+h-1, to be filled-in
//		int w,     // width of the image domain
//		int h,     // height of the image domain
//		float p,   // cosine of direction
//		float q    // sine of direction
//		)
//{
//	int N = w*h + 2*(w+h);
//	for (int i = 0; i < N; i++)
//		o[i] = -1;
//
//	if (fabs(p) <= q) // canonical case
//		canonical_bresenham_parkour(o, w, h, p, q);
//	else if (fabs(p) <= -q) { // swap the sign of q
//		canonical_bresenham_parkour(o, w, h, p, -q);
//		for (int i = 0; i < N; i++) if (o[i] >= 0)
//			o[i] = (h - 1 - o[i]/w)*w + o[i]%w;
//	} else if (fabs(q) < p) { // swap q and p
//		canonical_bresenham_parkour(o, h, w, q, p);
//		for (int i = 0; i < N; i++) if (o[i] >= 0)
//			o[i] = (o[i]%h)*w + o[i]/h;
//	} else if (fabs(q) < -p) { // swap q and -p
//		//fprintf(stderr, "here pq = %g %g!", p, q);
//		canonical_bresenham_parkour(o, h, w, q, -p);
//		for (int i = 0; i < N; i++) if (o[i] >= 0)
//			o[i] = (o[i]%h)*w + w - 1 - o[i]/h;
//	}
//}

#include "xmalloc.c"

static float getpix1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

static float bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static float getpixel_bilinear(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpix1(x, w, h, ip  , iq  );
	float b = getpix1(x, w, h, ip+1, iq  );
	float c = getpix1(x, w, h, ip  , iq+1);
	float d = getpix1(x, w, h, ip+1, iq+1);
	float r = bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

static void cast_shadows(
		float *D,      // DEM raster data, to be filled-in with NAN
		int w,         // width of the raster
		int h,         // height of the raster
		float p,       // cosine of ground-projected sun direction
		float q,       // sine of ground-projected sun direction
		float a        // slope of the sun direction
		)
{
	char *M = xmalloc(w*h*sizeof*M);
	for (int i = 0; i < w*h; i++) M[i] = 1;
	if (!p && !q) return;
	fprintf(stderr, "cast shadows pqa = %g %g %g\n", p, q, a);
	int N = w*h + 2*(w+h);  // maximum possible length of parkour
	int *i = xmalloc(N*sizeof*i); // indices of Bresenham parkour
	float *oxy = xmalloc(2*N*sizeof*oxy);
	compute_the_bresenham_parkour_hack(i, oxy, w, h, p, q);
	int l = -1; // index of the first point on the current line
	            // note: l is an index for the "i" and "oxy" arrays
	bool debbie = false;
	////float *dbuf = xmalloc(2*w*h*sizeof*dbuf);
	////for (int i = 0; i < 2*w*h; i++) dbuf[i] = 42;
	float xy0[2];
	for (int j = 0; j < N; j++)
	{
		if (i[j] < 0 && j+1<N && i[j+1]) // beginning of a line
		{
			l = j + 1;
			//if (i[l] > 0) debbie = i[l]%w==0 && i[l]/w==380;
			xy0[0] = i[l] % w;
			xy0[1] = i[l] / w;
		}
		if (i[j] < 0 || i[l] < 0) continue;
		assert(i[j] >= 0);
		assert(i[j] < w*h);
		if (i[l] < 0) fprintf(stderr, "SHIT l=%d i[l]=%d\n", l, i[l]);
		////dbuf[2*i[j]+0] = oxy[2*j+0];
		////dbuf[2*i[j]+1] = oxy[2*j+1];
		 int iX = i[l] % w;  // X of last occluding point
		 int iY = i[l] / w;  // Y of last occluding point
		 assert(0 <= iX && iX < w);
		 assert(0 <= iY && iY < h);
		float X = oxy[2*l+0];
		float Y = oxy[2*l+1];
		 if (fabs(X-iX)>0.6) fprintf(stderr, "iX=%d X=%g\n", iX, X);
		 assert(fabs(X-iX) < 0.6);
		 assert(fabs(Y-iY) < 0.6);
		//float Z = D[i[l]];   // Z of last occluding point
		float Z = getpixel_bilinear(D, w, h, X, Y);
		 int ix = i[j] % w;  // x of current point
		 int iy = i[j] / w;  // y of current point
		 assert(0 <= ix && ix < w);
		 assert(0 <= iy && iy < h);
		float x = oxy[2*j+0];
		float y = oxy[2*j+1];
		 assert(fabs(x-ix) < 0.6);
		 assert(fabs(y-iy) < 0.6);
		//float z = D[i[j]];   // z of current point
		float z = getpixel_bilinear(D, w, h, x, y);
		if (debbie)
			fprintf(stderr,
				"xy0 = %g %G  iXY=%d %d, XYZ=%g %g %g, "
				"ixy=%d %d, xyz = %g %g %g, " //xxyy=%g %g, "
				"j=%d, i[j]=%d\n",
				xy0[0], xy0[1], iX, iY, X, Y, Z,
				ix, iy, x, y, z, //xx, yy,
				j, i[j]);
		if (Z - z > a * hypot(x-X, y-Y)) // occluded point
		//if (Z - z > a * hypot(ix-iX, iy-iY)) // occluded point
			//D[i[j]] = NAN;
			M[i[j]] = 0;
		else // new occluding point
			l = j;
	}
	free(i);
	free(oxy);
	////iio_write_image_float_vec("/tmp/dbuf.npy", dbuf, w, h, 2);
	free(dbuf);
	for (int i = 0; i < w*h; i++)
		if (!M[i])
			D[i] = NAN;
	free(M);
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

#ifndef OMIT_MAIN_SHADOWCAST
//#define MAIN_VERTSHADOW
#define MAIN_THREEDSHADOW
#endif

#ifdef MAIN_VERTSHADOW
#include <stdio.h>
#include <stdlib.h>
//#include "iio.h"      // library for image input/output
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
//#include "iio.h"      // library for image input/output
#include "pickopt.c"  // function "pick_option" for processing args
int main(int c, char *v[])
{
	// process input arguments
	_Bool m = pick_option(&c, &v, "m", NULL);  // shadow mask
	_Bool M = pick_option(&c, &v, "M", NULL);  // negative shadow mask
	if (c < 4 || c > 6) {
		fprintf(stderr, "usage:\n\t%s p q r [dem_in [dem_out]]\n", *v);
		//                          0 1 2 3  4       5
		return 1;
	}
	float param_p = atof(v[1]); // p,q,r = direction of the sun
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

	// cast the shadows
	cast_shadows(x, w, h, param_p, param_q, param_r);

	// if mask is requested, create a binary mask
	if (m) for (int i = 0; i < w*h; i++)
		x[i] = 255*isnan(x[i]);
	else if (M) for (int i = 0; i < w*h; i++)
		x[i] = 255*!isnan(x[i]);

	// save the output image
	iio_write_image_float_split(filename_out, x, w, h, 1);

	// cleanup (unnecessary) and exit
	free(x);
	return 0;
}
#endif//MAIN_THREEDSHADOW


