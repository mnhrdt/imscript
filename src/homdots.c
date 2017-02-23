#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


static void apply_homography(double y[2], double H[9], double x[2])
{
	double X = H[0] * x[0] + H[1] * x[1] + H[2];
	double Y = H[3] * x[0] + H[4] * x[1] + H[5];
	double Z = H[6] * x[0] + H[7] * x[1] + H[8];
	y[0] = X / Z;
	y[1] = Y / Z;
}

static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

void homdots(float *out, int out_w, int out_h,
		float *x, int w, int h, int pd,
		double H[9])
{
	int *counter = malloc(out_w * out_h * sizeof*counter);
	for (int i = 0; i < out_w * out_h; i++)
		counter[i] = 0;

	for (int i = 0; i < out_w * out_h; i++)
	for (int l = 0; l < pd; l++)
		out[pd*i+l] = 0;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double p[2] = {i, j}, q[2];
		apply_homography(q, H, p);
		int qi = lrint(q[0]);
		int qj = lrint(q[1]);
		if (insideP(out_w, out_h, qi, qj))
		{
			int oidx = qj * out_w + qi;
			int xidx = j * w + i;
			counter[oidx] += 1;
			for (int l = 0; l < pd; l++)
				out[ pd * oidx + l ] += x[ pd * xidx + l ];
		}
	}

	for (int i = 0; i < out_w * out_h; i++)
	for (int l = 0; l < pd; l++)
		if (counter[i])
			out[pd*i+l] /= counter[i];
		else
			out[pd*i+l] = NAN;
}

#include <stdio.h>
#include "iio.h"

int main(int c, char *v[])
{
	// process input arguments
	if (c < 12 || c > 14 ) {
		fprintf(stderr, "usage:\n\t"
			"%s h1 ... h9 W H [in.png [out.tiff]]\n", *v);
		//       0  1      9  10 11 12      13
		return 1;
	}
	double H[9]; for (int i = 0; i < 9; i++) H[i] = atof(v[1+i]);
	int out_w = atoi(v[10]);
	int out_h = atoi(v[11]);
	char *filename_in  = c > 12 ? v[12] : "-";
	char *filename_out = c > 13 ? v[13] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *y = malloc(out_w * out_h * pd * sizeof*y);
	if (!y) return 2;

	// perform the computation
	homdots(y, out_w, out_h, x, w, h, pd, H);

	// save output result
	iio_write_image_float_vec(filename_out, y, out_w, out_h, pd);

	// cleanup and exit
	free(x);
	free(y);
	return 0;
}
