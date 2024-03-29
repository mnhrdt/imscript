// watermask: heuristic intended to find water planes in srtm4 data
// turns large connected components of constant value into NAN

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "iio.h"
#include "fail.c"
#include "xmalloc.c"
#include "ccproc.c"
#include <stdarg.h>

static void img_debug(float *x, int w, int h, int pd, const char *fmt, ...)
{
	va_list ap;
	char fname[FILENAME_MAX];
	va_start(ap, fmt);
	vsnprintf(fname, FILENAME_MAX, fmt, ap);
	va_end(ap);
	iio_write_image_float_vec(fname, x, w, h, pd);
}

static void img_debug_int(int *x, int w, int h, int pd, const char *fmt, ...)
{
	va_list ap;
	char fname[FILENAME_MAX];
	va_start(ap, fmt);
	vsnprintf(fname, FILENAME_MAX, fmt, ap);
	va_end(ap);
	iio_write_image_int_vec(fname, x, w, h, pd);
}
static bool size_is_admissible(double *size, int nsizes, int sidx, int x)
{
	assert(size[nsizes-1] == INFINITY);
	if (sidx == 0)
		return x <= size[sidx];
	if (sidx == nsizes-1)
		return size[sidx-1] < x;
	return size[sidx-1] < x && x <= size[sidx];
}

int main(int c, char *v[])
{
	if (c != 2 && c  != 3 && c != 4)
		return fprintf(stderr, "usage:\n\t%s n [in [out]]\n", *v);
	//                                         0 1  2   3
	int minsize = atoi(v[1]);
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (pd != 1) return fprintf(stderr, "only gray images please\n");
	float *y = xmalloc(w*h*sizeof*y);

	int *out_size   = xmalloc(w*h*sizeof*out_size);
	int *out_bdsize = xmalloc(w*h*sizeof*out_size);
	int *out_all    = xmalloc(w*h*sizeof*out_size);
	int *out_first  = xmalloc(w*h*sizeof*out_size);
	int *out_idx    = xmalloc(w*h*sizeof*out_size);

	int r = ccproc(out_size, out_bdsize, out_all, out_first, out_idx,
			x, w, h,
			//float_eq_isnan
			floatnan_equality
			//floats_are_equal
			);

	fprintf(stderr, "r = %d\n", r);

	for (int i = 0; i < w*h; i++)
		y[i] = out_size[ out_idx[i] ] < minsize ? x[i] : NAN;

	iio_write_image_float(filename_out, y, w, h);

	return 0;
	//for (int i = 0; i < r; i++)
	//{
	//	fprintf(stderr, "r_%d : area=%d perim=%d first=%d\n",
	//			i, out_size[i], out_bdsize[i], out_first[i]);
	//	fprintf(stderr, "\tbd :");
	//	for (int j = 0; j < out_bdsize[i]; j++)
	//		fprintf(stderr, " %d", out_all[out_first[i]+j]);
	//	fprintf(stderr, "\n");
	//	fprintf(stderr, "\tin :");
	//	for (int j = 0; j < out_size[i] - out_bdsize[i]; j++)
	//		fprintf(stderr, " %d", out_all[out_first[i]+out_bdsize[i]+j]);
	//	fprintf(stderr, "\n");
	//}

	////for (int s = 0; s < nsizes; s++)
	////{
	////	int counts = 0;
	////	for (int i = 0; i < w*h; i++)
	////		x[i] = -1000;
	////	for (int i = 0; i < r; i++)
	////		if (size_is_admissible(size, nsizes, s, out_size[i]))
	////		{
	////			for (int j = 0; j < out_size[i]; j++)
	////				x[out_all[out_first[i]+j]] = i;
	////			counts += 1;
	////		}
	////	fprintf(stderr, "count upto %g = %d\n", size[s], counts);
	////	img_debug(x, w, h, 1, "cosa_upto_%g.tiff", size[s]);
	////}

	//fprintf(stderr, "ccproc returned %d\n", r);
	//for (int i = 0; i < r; i++)
	//{
	//	fprintf(stderr, "out_size[%d/%d] = %d (%d)\n",
	//			i, r, out_size[i], out_bdsize[i]);
	//	if (i > 100) break;
	//}

	//// verify consistence
	//int totsize = 0;
	//for (int i = 0; i < r; i++)
	//{
	//	for (int j = 0; j < out_size[i]; j++)
	//		x[out_all[out_first[i] + j]] = i;
	//	totsize += out_size[i];
	//}
	//assert(totsize == w*h);

	//iio_write_image_int("ccproc_idx.tiff", out_idx, w, h);
	//iio_write_image_int("ccproc_all.tiff", out_all, w, h);
	//iio_write_image_float("ccproc_xxx.tiff", x, w, h);

	//for (int i = 0; i < w*h; i++)
	//	x[i] = -1;
	//for (int i = 0; i < r; i++)
	//{
	//	for (int j = 0; j < out_bdsize[i]; j++)
	//		x[out_all[out_first[i] + j]] = i;
	//}
	//iio_write_image_float("ccproc_yyy.tiff", x, w, h);


	//free(out_size);
	//free(out_bdsize);
	//free(out_all);
	//free(out_first);
	//free(out_idx);
	//free(x);
	return 0;
}
