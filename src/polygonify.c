#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include "ccproc.c"

#include "smapa.h"
SMART_PARAMETER(POLYGONIFY_a,50)
SMART_PARAMETER(POLYGONIFY_A,30000)
SMART_PARAMETER(POLYGONIFY_F,2)


static double transport_distance(double *a, double *b, int n)
{
	// accumulate
	double aa[n], bb[n];
	aa[0] = a[0];
	bb[0] = b[0];
	for (int i = 1; i < n; i++)
	{
		aa[i] = aa[i-1] + a[i];
		bb[i] = bb[i-1] + b[i];
	}

	// normalize
	for (int i = 0; i < n; i++)
	{
		aa[i] /= aa[n-1];
		bb[i] /= bb[n-1];
	}

	// compute L1 norm of the difference
	double r = 0;
	for (int i = 0; i < n; i++)
		r += fabs(aa[i] - bb[i]);

	return r;
}


#include "iio.h"
int main(int c, char *v[])
{
	if (c != 5)
		return fprintf(stderr, "usage:\n\t"
		"%s M1.tif M2.tif 11band.tif segm.tif > polys.txt\n", *v);
	//        0 1      2      3          4
	char *filename_M1 = v[1];
	char *filename_M2 = v[2];
	char *filename_x = v[3];
	char *filename_y = v[4];

	int wm1, hm1, wm2, hm2, wx, hx, pd, wy, hy;
	double *M1 = iio_read_image_double(filename_M1, &wm1, &hm1); //histo_in
	double *M2 = iio_read_image_double(filename_M2, &wm2, &hm2); //histo_out
	float *x = iio_read_image_float_vec(filename_x, &wx, &hx, &pd); //11band
	float *y = iio_read_image_float(filename_y, &wx, &hx); // segmentation
	assert(wm1 == 256);
	assert(wm2 == 256);
	assert(hm1 == 11);
	assert(hm2 == 11);
	assert(pd == 11);
	assert(wx == wy);
	assert(hx == hy);
	int w = wx;
	int h = hx;
	int nbins = 256;
	int nbands = 11;

	int *out_size   = xmalloc(w*h*sizeof*out_size);
	int *out_bdsize = xmalloc(w*h*sizeof*out_size);
	int *out_all    = xmalloc(w*h*sizeof*out_size);
	int *out_first  = xmalloc(w*h*sizeof*out_size);
	int *out_idx    = xmalloc(w*h*sizeof*out_size);
	int *bd_tmp     = xmalloc(3*4*w*h*sizeof*bd_tmp);

	int r = ccproc(out_size, out_bdsize, out_all, out_first, out_idx,
			y, w, h, floatnan_equality);

	for (int i = 0; i < r; i++)
	{
		// reject region if it is too small or too big
		if (out_size[i] < POLYGONIFY_a())
			continue;
		if (out_size[i] > POLYGONIFY_A())
			continue;

		// compute 11-hitogram inside this region
		double H[nbins*nbands];
		for (int l = 0; l < nbins*nbands; l++)
			H[l] = 0;
		for (int j = 0; j < out_size[i]; j++)
		{
			int pidx = out_all[out_first[i]+j];
			float *v = x + nbands * pidx;
			for (int l = 0; l < nbands; l++)
			{
				float vl = v[l];
				if (!isfinite(vl)) continue;
				int bin = floor(vl);
				if (bin < 0) bin = 0;
				if (bin >= nbins) bin = nbins - 1;
				H[l*nbins+bin] += 1;
			}
		}

		// compute distances to both barycenters
		double d1 = transport_distance(H, M1, nbins*nbands);
		double d2 = transport_distance(H, M2, nbins*nbands);
		fprintf(stderr, "region %d\tarea %d\td1=%g\td2=%g\tr=%g\n", i,
				out_size[i], d1, d2, d2/d1);

		// if it is much closer to M1 than to M2, output the polygon
		if (POLYGONIFY_F() * d1 < d2)
		{
			fprintf(stderr, "\taccepted!\n");
			// index of one interior pixel
			int ridx = out_all[out_first[i]+0];
			int rx = ridx % w;
			int ry = ridx / w;
			int bd_n = bfollow(bd_tmp, y, w, h, floatnan_equality,
					rx, ry);
			for (int k = 0; k < bd_n; k++)
			{
				float px = bd_tmp[3*k+0];
				float py = bd_tmp[3*k+1];
				switch(bd_tmp[3*k+2]) {
				case 0: px += 0.5; break;
				case 1: py -= 0.5; break;
				case 2: px -= 0.5; break;
				case 3: py += 0.5; break;
				}
				printf("%g %g ", px, py);
			}
			printf("%g\n", d2/d1);
		}
	}
	return 0;
}
