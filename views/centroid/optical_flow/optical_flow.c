#include <math.h>
#include <stdio.h>
#include <stdlib.h>



// this function is like "malloc", but it returns always a valid pointer
static void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		exit(fprintf(stderr, "out of memory\n"));
	return p;
}



// compute a smoothness measure of a 3x3 patch from a vector field (U,V)
//   [0   , 0.5] = very smooth
//   [0.5 ,   1] = smooth enough
//   [1   , ...] = not smooth
static float patch_smoothness(float U[3][3], float V[3][3])
{
	// upwind directional derivatives
	float Ux = fmax(fabs(U[1][2] - U[1][1]), fabs(U[1][1] - U[1][0]));
	float Uy = fmax(fabs(U[2][1] - U[1][1]), fabs(U[1][1] - U[0][1]));
	float Vx = fmax(fabs(V[1][2] - V[1][1]), fabs(V[1][1] - V[1][0]));
	float Vy = fmax(fabs(V[2][1] - V[1][1]), fabs(V[1][1] - V[0][1]));

	// laplacian
	float lU = fabs(4*U[1][1] - U[2][1] - U[1][2] - U[0][1] - U[1][0]);
	float lV = fabs(4*V[1][1] - V[2][1] - V[1][2] - V[0][1] - V[1][0]);

	// divergence
	float Uxx = -2 * U[1][1] + U[1][2] + U[1][0];
	float Vyy = -2 * V[1][1] + V[2][1] + V[0][1];
	float dUV = fabs(Uxx + Vyy);

	// max of all the above
	return fmax(fmax(fmax(fmax(fmax(fmax(Ux,Uy),Vx),Vy),lU),lV),dUV);
}

typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by nearest value
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// s := pointwise smoothness of the vector field (u,v)
static void field_smoothness(float *s, float *u, float *v, int w, int h)
{
	getpixel_operator p = getpixel_1;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float U[3][3], V[3][3];
		for (int jj = 0; jj < 3; jj++)
		for (int ii = 0; ii < 3; ii++)
		{
			U[jj][ii] = p(u, w, h, i+ii-1, j+jj-1);
			V[jj][ii] = p(v, w, h, i+ii-1, j+jj-1);
		}
		s[w*j+i] = patch_smoothness(U, V);;
	}
}

#define SMOOTHNESS_THRESHOLD 1.1
#define BAD_NEIGHBORHOOD_RADIUS 1

// mark the pixel (i,j) and its neighbors as bad (indicated by NAN)
static void mark_bad_neighborhood(float *x, int w, int h, int i, int j)
{
	int radius = BAD_NEIGHBORHOOD_RADIUS;
	for (int dx = -radius; dx <= radius; dx++)
	for (int dy = -radius; dy <= radius; dy++)
	{
		int ii = i + dx;
		int jj = j + dy;
		if (ii >= 0 && ii < w && jj >= 0 && jj < h)
			x[jj*w + ii] = NAN;
	}
}

// (U,V) := (u,v) with smoothened discontinuities
static void field_smoothen(float *U, float *V, float *u, float *v, int w, int h)
{
	// find and mark bad pixels
	float *smoothness = xmalloc(w * h * sizeof(float));
	field_smoothness(smoothness, u, v, w, h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (fabs(smoothness[w*j+i]) > SMOOTHNESS_THRESHOLD)
		{
			mark_bad_neighborhood(u, w, h, i, j);
			mark_bad_neighborhood(v, w, h, i, j);
		}
	free(smoothness);

	// count number of masked pixels (only for reporting)
	int nnan = 0;
	for (int i = 0; i < w*h; i++)
		if (!isfinite(u[i]))
			nnan += 1;
	fprintf(stderr, "\tmasked=%g%%", nnan*100.0/(w*h));

	// fill-in the holes
	float timestep = 0.25;
	int niter=10, nscales=10;
	void elap_recursive(float*,float*,int,int,float,int,int);
	elap_recursive(U, u, w, h, timestep, niter, nscales);
	elap_recursive(V, v, w, h, timestep, niter, nscales);
}

// convert a multispectral image to gray
static void get_gray_values(float *gx, float *x, int w, int h, int pd)
{
	for (int i = 0; i < w*h; i++)
	{
		gx[i] = 0;
		for (int j = 0; j < pd; j++)
			gx[i] += x[i*pd+j]/pd;
	}
}

// u := optical flow from image a to image b
void optical_flow(float *u, float *a, float *b, int w, int h, int pd)
{
	// temporary stockage
	float *gray_a = xmalloc(w * h * sizeof(float));
	float *gray_b = xmalloc(w * h * sizeof(float));
	float *flow_x = xmalloc(2 * w * h * sizeof(float));
	float *flow_y = flow_x + w*h;

	// prepare input data
	get_gray_values(gray_a, a, w, h, pd);
	get_gray_values(gray_b, b, w, h, pd);

	// set Horn-Schunck parameters
	float alpha = 20;
	int nscales = 10;
	float zfactor = 0.5;
	int nwarps = 2;
	float epsilon = 0.1;
	int maxiter = 50;
	int verbose = 0;

	// compute optical flow
	void horn_schunck_pyramidal(float*,float*,float*,float*,int,int,
			float,int,float,int,float,int,int);
	horn_schunck_pyramidal(gray_a, gray_b, flow_x, flow_y, w, h,
			alpha, nscales, zfactor,
			nwarps, epsilon, maxiter,
			verbose);

	// smoothen flow
	float *smooth_x = xmalloc(2 * w * h * sizeof(float));
	float *smooth_y = smooth_x + w*h;
	field_smoothen(smooth_x, smooth_y, flow_x, flow_y,  w, h);

	// prepare output data
	for (int i = 0; i < w*h; i++)
	{
		u[2*i+0] = smooth_x[i];
		u[2*i+1] = smooth_y[i];
	}

	// cleanup
	free(smooth_x);
	free(gray_a);
	free(gray_b);
	free(flow_x);
}
