#define MAIN_CENTROID


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// small inline function for bicubic interpolation
#include "bicubic.c"


// naming convention for all the functions below:
// w: width of the image domain
// h: height of the image domain
// pd: pixel dimension (e.g., 1 or 3)


// y := pull back of image x by vector field f
static void field_pull_back(float *y, float *f, float *x, int w, int h, int pd)
{
	float (*flow)[w][2] = (void*)f;
	float (*out)[w][pd] = (void*)y;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p[2] = {i + flow[j][i][0], j + flow[j][i][1]};
		bicubic_interpolation(out[j][i], x, w, h, pd, p[0], p[1]);
	}
}

// v(x) := -u(x+v(x))
static void inversion_iteration(float *v, float *u, int w, int h)
{
	float (*V)[w][2] = (void*)v;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p[2];
		float qx = i + V[j][i][0];
		float qy = j + V[j][i][1];
		bicubic_interpolation(p, u, w, h, 2, qx, qy);
		V[j][i][0] = -p[0];
		V[j][i][1] = -p[1];
	}
}

// v := inverse vector field of u
static void field_invert(float *v, float *u, int w, int h)
{
	for (int i = 0; i < w*h*2; i++)
		v[i] = 0;

	int niter = 7;
	for (int i = 0; i < niter; i++)
		inversion_iteration(v, u, w, h);

}


// f := optical flow between images a and b
// (this function shall be defined elsewhere)
void optical_flow(float *f, float *a, float *b, int w, int h, int pd);

// auxiliary function (like malloc, but always returns a valid pointer)
static void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		exit(fprintf(stderr,
			"ERROR: out of mem requesting %zu butes\n", size));
	return p;
}

// out := centroid of images in[n] from image ref
void centroid(float *out, float **in, int n, float *ref, int w, int h, int pd)
{
	// storage for optical flows
	float *flow[n];
	for (int i = 0; i < n; i++)
		flow[i] = xmalloc(w * h * 2 * sizeof(float));
	float *flow_avg = xmalloc(w * h * 2 * sizeof(float));
	float *flow_iavg= xmalloc(w * h * 2 * sizeof(float));

	// compute all optical flows
	for (int i = 0; i < n; i++)
	{
		fprintf(stderr, "flow %d/%d", i+1, n);
		optical_flow(flow[i], ref, in[i], w, h, pd);
		fprintf(stderr, "\n");
	}

	// compute average
	for (int i = 0; i < 2*w*h; i++)
		flow_avg[i] = 0;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < 2*w*h; i++)
			flow_avg[i] += flow[j][i] / n;

	// push-forward the reference image by average flow
	fprintf(stderr, "warping...\n");
	field_invert(flow_iavg, flow_avg, w, h);
	field_pull_back(out, flow_iavg, ref, w, h, pd);

	// cleanup
	free(flow_avg);
	free(flow_iavg);
	for (int i = 0; i < n; i++)
		free(flow[i]);
}


#ifdef MAIN_CENTROID
#include "iio.h"

int main(int c, char *v[])
{
	// process input arguments
	if (c != 6) {
		fprintf(stderr, "usage:\n\t"
				"%s in_pattern first last refimg out\n", *v);
		//                0 1          2     3    4      5
		return 1;
	}
	char *filename_in_fmt = v[1];
	int index_first = atoi(v[2]);
	int index_last = atoi(v[3]);
	char *filename_ref = v[4];
	char *filename_out = v[5];

	// compute number of input frames
	int n = index_last - index_first + 1;
	if (n < 2)
		exit(fprintf(stderr,"need at least 2 images (got %d)\n", n));

	// read input images
	int w, h, pd;
	float *ref, *in[n];
	ref = iio_read_image_float_vec(filename_ref, &w, &h, &pd);
	for (int i = 0; i < n; i++)
	{
		char fname[FILENAME_MAX];
		snprintf(fname, FILENAME_MAX, filename_in_fmt, i + index_first);
		int ww, hh, ppdd;
		in[i] = iio_read_image_float_vec(fname, &ww, &hh, &ppdd);
		if (w != ww || h != hh || pd != ppdd)
			exit(fprintf(stderr,"input sizes mismatch "
					"(%d,%d,%d)!=(%d,%d,%d)",
					w, h, pd, ww, hh, ppdd));
	}

	// allocate space for output image
	float *out = xmalloc(w * h * pd * sizeof*out);

	// compute centroid
	centroid(out, in, n, ref, w, h, pd);

	// save output image
	iio_save_image_float_vec(filename_out, out, w, h, pd);

	// free resources and exit
	for (int i = 0; i < n; i++)
		free(in[i]);
	free(ref);
	free(out);
	return 0;
}
#endif//MAIN_CENTROID
