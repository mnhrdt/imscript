#include <math.h>
#include <stdlib.h>

static void *xmalloc(size_t n)
{
	void *p = malloc(n);
	if (!p)
		abort();
	return p;
}

static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

static float laplacian(float *x, int w, int h, int i, int j)
{
	return  -4 * getpixel(x, w, h, i  , j  )
		     + getpixel(x, w, h, i+1, j  )
		     + getpixel(x, w, h, i  , j+1)
		     + getpixel(x, w, h, i-1, j  )
		     + getpixel(x, w, h, i  , j-1);
}

static void perform_one_iteration(float *x, float *dat,
		int *mask, int nmask, int w, int h, float tstep)
{
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[2*p+0];
		int j = mask[2*p+1];
		int idx = j*w + i;
		x[idx] += tstep * (laplacian(x, w, h, i, j) - dat[idx]);
	}
}

static int *build_mask(int *out_nmask, float *x, int w, int h)
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int *mask = xmalloc(w*h*2*sizeof*mask), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i]))
		{
			mask[2*cx+0] = i;
			mask[2*cx+1] = j;
			cx += 1;
		}
	*out_nmask = nmask;
	return mask;
}

static void poisson_extension_with_init(float *out, float *inb, float *dat,
		int w, int h,
		float timestep, int niter, float *initialization)
{
	int nmask, *mask = build_mask(&nmask, inb, w, h);
	for (int i = 0; i < w*h; i++)
		out[i] = isfinite(inb[i]) ? inb[i] : initialization[i];
	for (int i = 0; i < niter; i++)
		perform_one_iteration(out, dat, mask, nmask, w, h, timestep);
	free(mask);
}

static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = getpixel(in, iw, ih, 2*i, 2*j);
		a[1] = getpixel(in, iw, ih, 2*i+1, 2*j);
		a[2] = getpixel(in, iw, ih, 2*i, 2*j+1);
		a[3] = getpixel(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k]))
			{
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
}

static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float x = round((i - 0.5)/2);
		float y = round((j - 0.5)/2);
		out[ow*j+i] = getpixel(in, iw, ih, x, y);
	}
}

void poisson_rec(float *out, float *in, float *dat, int w, int h,
		float tstep, int niter, int scale)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins  = xmalloc(ws * hs * sizeof*ins);
		float *dats = xmalloc(ws * hs * sizeof*dats);
		float *outs = xmalloc(ws * hs * sizeof*outs);
		zoom_out_by_factor_two(ins, ws, hs, in, w, h);
		zoom_out_by_factor_two(dats, ws, hs, dat, w, h);
		poisson_rec(outs, ins, dats, ws, hs, tstep, niter, scale-1);
		zoom_in_by_factor_two(init, w, h, outs, ws, hs);
		free(ins);
		free(dats);
		free(outs);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	poisson_extension_with_init(out, in, dat, w, h, tstep, niter, init);
	free(init);
}
