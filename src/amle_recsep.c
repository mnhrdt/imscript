#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// utility function that always returns a valid pointer to memory
static void *xmalloc(size_t n)
{
	void *new = malloc(n);
	if (!new)
	{
		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
		exit(1);
	}
	return new;
}


// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
inline static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w];
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
inline static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}


// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[2]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
		}
	assert(cx == nmask);

	*out_nmask = nmask;
	return mask;
}

#include "smapa.h"
SMART_PARAMETER_SILENT(AMLE_NN,4)

// evaluate an image in a neighborhood
static int get_nvals(float *v, float *wv2, float *x, int w, int h, int i, int j)
{
	int r = 0, n[][3] = { // {x, y, x*x+y*y}
		{+1,0,1}, {0,+1,1}, {-1,0,1}, {0,-1,1}, // 4-connexity
		{+1,+1,2}, {-1,-1,2}, // 6-connexity
		{-1,+1,2}, {+1,-1,2}, // 8-connexity
		{+2,+1,5}, {+1,+2,5}, {+2,-1,5}, {+1,-2,5},
		{-2,-1,5}, {-1,-2,5}, {-2,+1,5}, {-1,+2,5}, // 16-neighbors
		{+3,+1,10}, {+1,+3,10}, {+3,-1,10}, {+1,-3,10},
		{-3,-1,10}, {-1,-3,10}, {-3,+1,10}, {-1,+3,10}, // 24-neighbors
		{+3,+2,13}, {+2,+3,13}, {+3,-2,13}, {+2,-3,13},
		{-3,-2,13}, {-2,-3,13}, {-3,+2,13}, {-2,+3,13}, // 32-neighbors
		//{+2,0,4}, {0,+2,4}, {-2,0,4}, {0,-2,4}, (non-primitive)
		//{+4,0,16}, {0,+4,16}, {-4,0,16}, {0,-4,16} (non-primitive)
	};
	int n33[32][3] = {
		{+1,0,1}, {0,+1,1}, {-1,0,1}, {0,-1,1}, // 4-connexity
		{+1,+1,2}, {-1,-1,2}, // 6-connexity
		{-1,+1,2}, {+1,-1,2}, // 8-connexity
		{+4,0,16}, {-4,0,16}, {0,4,16}, {0,-4,16},//4
		{+4,+1,17}, {-4,+1,17}, {+4,-1,17}, {-4,-1,17},
		{+1,+4,17}, {-1,+4,17}, {+1,-4,17}, {-1,-4,17},//12
		{+4,+2,20}, {-4,+2,20}, {+4,-2,20}, {-4,-2,20},
		{+2,+4,20}, {-2,+4,20}, {+2,-4,20}, {-2,-4,20},
		{+3,+3,18}, {-3,+3,18}, {+3,-3,18}, {-3,-3,18},//24
		//{+4,+3,25}, {-4,+3,25}, {+4,-3,25}, {-4,-3,25},
		//{+3,+4,25}, {-3,+4,25}, {+3,-4,25}, {-3,-4,25},//28
		//{+4,+4,32}, {-4,+4,32}, {+4,-4,32}, {-4,-4,32},//32

		//{+3,+2,13}, {-3,+2,13}, {+3,-2,13}, {-3,-2,13},
		//{+2,+3,13}, {-2,+3,13}, {+2,-3,13}, {-2,-3,13},
		//{+4,+3,25}, {-4,+3,25}, {+4,-3,25}, {-4,-3,25},
		//{+3,+4,25}, {-3,+4,25}, {+3,-4,25}, {-3,-4,25},
	};
	int (*pn)[3] = n;
	int nn = AMLE_NN();
	if (nn == 22) { pn = n33; nn = 32; }
	if (nn == 33) { pn = n33; nn = 32; }
	if (nn == 44) { pn = n33; nn = 40; }
	for (int p = 0; p < nn; p++)
	{
		int ii = i + pn[p][0];
		int jj = j + pn[p][1];
		if (ii >= 0 && jj >= 0 && ii < w && jj < h)
		{
			v[r] = x[w*jj+ii];
			if (wv2)
				//wv2[r] = pow(hypot(pn[p][0], pn[p][1]),2);
				wv2[r] = pn[p][2];
			r += 1;
		}
	}
	return r;
}

static void get_minmax(float *min, float *max, float *x, int n)
{
	*min = INFINITY;
	*max = -INFINITY;
	for (int i = 0; i < n; i++)
	{
		if (x[i] < *min) *min = x[i];
		if (x[i] > *max) *max = x[i];
	}
}

static void get_minmax_idx(int *min, int *max, float *x, int n)
{
	*min = *max = 0;
	for (int i = 1; i < n; i++)
	{
		if (x[i] < x[*min]) *min = i;
		if (x[i] > x[*max]) *max = i;
	}
}

static float amle_iteration(float *x, int w, int h, int (*mask)[2], int nmask)
{
	float actus = 0;
	float actumax = 0;
//#pragma omp parallel for
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i, min, max;
		float value[0x100], weight[0x100];
		int nv = get_nvals(value, weight, x, w, h, i, j);
		get_minmax_idx(&min, &max, value, nv);
		float a = sqrt(weight[max]);
		float b = sqrt(weight[min]);
		// TODO: clarify the appropriate normalization
		//float a = weight[max];
		//float b = weight[min];
		float newx = (a*value[min] + b*value[max]) / (a + b);
		actus += fabs(x[idx] - newx);
		if (fabs(x[idx]-newx) > actumax)
			actumax = fabs(x[idx]-newx);
		x[idx] = newx;
	}
	return actumax;
}

// fill the holes of the image x using an infinity harmonic function
static void inf_harmonic_extension_with_init(
		float *y,        // output image
		float *x,        // input image (NAN values indicate holes)
		int w,           // image width
		int h,           // image height
		int niter,       // number of iterations to run
		float *initialization
		)
{
	// build list of masked pixels
	int nmask, (*mask)[2] = build_mask(&nmask, x, w, h);

	// initialize the solution to the given data at the masked pixels
	for (int i = 0; i < w*h; i++)
		y[i] = isfinite(x[i]) ? x[i] : initialization[i];

	// do the requested iterations
	for (int i = 0; i < niter; i++)
	{
		float u = amle_iteration(y, w, h, mask, nmask);

		//if (0 == i % 10)
		//fprintf(stderr, "size = %dx%d, i = %d, u = %g\n", w, h, i, u);
	}

	free(mask);
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i, 2*j);
		a[1] = p(in, iw, ih, 2*i+1, 2*j);
		a[2] = p(in, iw, ih, 2*i, 2*j+1);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

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

static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

static float bicubic_interpolation(float *img, int w, int h, float x, float y)
{
	x -= 1;
	y -= 1;

	getpixel_operator p = getpixel_1;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, ix + i, iy + j);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}



// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = p(in, iw, ih, round((i-0.5)/2), round((j-0.5)/2));
		//out[ow*j+i] = bicubic_interpolation(in, iw, ih, (i-0.5)/2.0, (j-0.5)/2.0);
		//out[ow*j+i] = bilinear_interpolation(in, iw, ih, (i-0.5)/2.0, (j-0.5)/2.0);
}

#include "iio.h"

SMART_PARAMETER_SILENT(AMLE_ONLY,0)

void amle_recursive(float *out, float *in, int w, int h, int niter, int scale)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins = xmalloc(ws*hs*sizeof*ins);
		float *outs = xmalloc(ws*hs*sizeof*outs);
		zoom_out_by_factor_two(ins, ws, hs, in, w, h);
		amle_recursive(outs, ins, ws, hs, niter, scale - 1);
		zoom_in_by_factor_two(init, w, h, outs, ws, hs);

		//char buf[FILENAME_MAX];
		//snprintf(buf, FILENAME_MAX, "/tmp/elap_rec_in_%d", scale);
		//iio_save_image_float(buf, in, w, h);
		//snprintf(buf, FILENAME_MAX, "/tmp/elap_rec_ins_%d", scale);
		//iio_save_image_float(buf, ins, ws, hs);
		//snprintf(buf, FILENAME_MAX, "/tmp/elap_rec_outs_%d", scale);
		//iio_save_image_float(buf, outs, ws, hs);
		//snprintf(buf, FILENAME_MAX, "/tmp/elap_rec_init_%d", scale);
		//iio_save_image_float(buf, init, w, h);

		free(ins);
		free(outs);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	if (AMLE_ONLY() > 0 && AMLE_ONLY()!=w) niter = 0;
	inf_harmonic_extension_with_init(out, in, w, h, niter, init);
	free(init);
}

// extension by AMLE each channel of a color image
void amle_recursive_separable(float *out, float *in, int w, int h, int pd,
		int niter, int nscales)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		amle_recursive(outl, inl, w, h, niter, nscales);
	}
}


int main(int argc, char *argv[])
{
	if (argc != 6) {
		fprintf(stderr, "usage:\n\t"
		"%s NITER NS data.png mask.png out.png\n", *argv);
		//0 1     2  3        4        5
		return 1;
	}
	int niter = atoi(argv[1]);
	int nscales = atoi(argv[2]);
	char *filename_in = argv[3];
	char *filename_mask = argv[4];
	char *filename_out = argv[5];

	int w[2], h[2], pd;
	float *in = iio_read_image_float_split(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch");
	float *out = xmalloc(*w**h*pd*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			for (int l = 0; l < pd; l++)
				in[*w**h*l+i] = NAN;

	amle_recursive_separable(out, in, *w, *h, pd, niter, nscales);

	iio_save_image_float_split(filename_out, out, *w, *h, pd);

	return 0;
}
