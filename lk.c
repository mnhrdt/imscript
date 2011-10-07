#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

static float sqr(float x)
{
	return x * x;
}

static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "out of memory\n"));
	return new;
}

typedef float (*extension_operator_float)(float*, int, int, int, int);

static float extend_float_image_by_zero(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i > w-1 || j > h-1)
		return 0;
	else
		return x[j*w+i];
}

static float extend_float_image_constant(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static void setpixel_float_image_vec(float *x, int w, int h, int pd,
		int i, int j, int l, float v)
{
	if (i < 0 || j < 0 || l < 0 || i>=w || j>=h || l >= pd)
		exit(fprintf(stderr, "bad setpixel (%d %d %d)[%d %d %d]\n",
					w, h, pd, i, j, l));
	float (*xx)[w][pd] = (void*)x;
	xx[j][i][l] = v;
}

static void compute_input_derivatives(float *Ex, float *Ey, float *Et,
		float *a, float *b, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		Ey[j*w+i] = (1.0/4) * ( p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j));
		Ex[j*w+i] = (1.0/4) * ( p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1));
		Et[j*w+i] = (1.0/4) * ( p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1));
	}
}

static void compute_bar(float *ubar, float *u, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ubar[j*w+i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}

static void hs_iteration(float *u, float *v,
		float *Ex, float *Ey, float *Et, int w, int h, float alpha)
{
	float *ubar = xmalloc(w * h * sizeof(float));
	float *vbar = xmalloc(w * h * sizeof(float));
	compute_bar(ubar, u, w, h);
	compute_bar(vbar, v, w, h);
	for (int i = 0; i < w*h; i++) {
		float t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		u[i] = ubar[i] - Ex[i] * t;
		v[i] = vbar[i] - Ey[i] * t;
	}
	free(ubar);
	free(vbar);
}

static void hs(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (int i = 0; i < w*h; i++)
		u[i] = v[i] = 0;
	for (int i = 0; i < niter; i++)
		hs_iteration(u, v, gx, gy, gt, w, h, alpha);
	free(gx);
	free(gy);
       	free(gt);
}

static void fill_window_offsets(int (*wo)[2], int kside)
{
	int idx = 0;
	int kradius = (kside - 1)/2;
	for (int j = 0; j < kside; j++)
	for (int i = 0; i < kside; i++)
	{
		wo[idx][0] = i - kradius;
		wo[idx][1] = j - kradius;
		idx += 1;
	}
}

static void fill_window_values(float *wv, int (*wo)[2], int kside, float sigma)
{
	int n = kside*kside;
	for (int i = 0; i < n; i++)
	{
		if (sigma < 0)
			wv[i] = 1;
		else {
			float r = hypot(wo[i][0], wo[i][1]);
			wv[i] = exp(-sqr(r/sigma));
		}
	}
	float m = 0;
	for (int i = 0; i < n; i++)
		m += wv[i];
	for (int i = 0; i < n; i++)
		wv[i] /= m;

	if (1) {
		FILE *f = fopen("/tmp/lk.kkk", "w");
		fprintf(f, "P2\n%d %d\n65535\n", kside, kside);
		for (int i = 0; i < n; i++)
			fprintf(f, "%g\n", wv[i]);
		fclose(f);
	}
}

//#include "svd.c"
#include "vvector.h"

static float solve_sdp_2x2(float x[2], float A[3], float b[2])
{
	float m[2][2] = {{A[0], A[1]}, {A[1], A[2]}};
	float n[2][2], det;
	INVERT_2X2(n, det, m);
	x[0] = n[0][0] * b[0] + n[0][1] * b[1];
	x[1] = n[1][0] * b[0] + n[1][1] * b[1];
	if(!det) {x[0]=x[1]=0;}
	//fprintf(stderr, "A=(%g %g %g), b=(%g %g), x=(%g %g)\n",
	//		A[0], A[1], A[2], b[0], b[1], x[0], x[1]);
	//float e[2];
	//e[0] = A[0]*x[0] + A[1]*x[1] - b[0];
	//e[1] = A[1]*x[0] + A[2]*x[1] - b[1];
	//fprintf(stderr, "e=(%g %g)\n", e[0], e[1]);
	return det;
}

#define STLEN 3

static void compute_structure_tensor_here(float atwa[STLEN],
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, int w, int h,
		int i, int j)
{
	extension_operator_float p = extend_float_image_constant;

	int n = kside * kside;
	atwa[0] = atwa[1] = atwa[2] = 0;

	for (int k = 0; k < n; k++)
	{
		int ii = i + wo[k][0];
		int jj = j + wo[k][1];
		atwa[0] += wv[k] * sqr(p(gx, w, h, ii, jj));
		atwa[1] += wv[k] * p(gx,w,h, ii, jj) * p(gy,w,h, ii, jj);
		atwa[2] += wv[k] * sqr(p(gy, w, h, ii, jj));
	}

	// svd of structure tensor (for visualization)
	//float A[4] = {atwa[0], atwa[1], atwa[1], atwa[2]};
	//float U[4], V[4], D[2];
	//svd_float(D, A, U, 2, V, 2);
	//atwa[3] = D[0];
	//atwa[4] = D[1];
	//atwa[5] = U[0];
	//atwa[6] = U[1];
	//atwa[7] = U[2];
	//atwa[8] = U[3];
}

static void compute_rhs_here(float rhsh[2],
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, float *gt, int w, int h,
		int i, int j)
{
	extension_operator_float p = extend_float_image_constant;

	int n = kside * kside;
	rhsh[0] = rhsh[1] = 0;

	for (int k = 0; k < n; k++)
	{
		int ii = i + wo[k][0];
		int jj = j + wo[k][1];
		rhsh[0] -= wv[k] * p(gx,w,h, ii, jj) * p(gt,w,h, ii, jj);
		rhsh[1] -= wv[k] * p(gy,w,h, ii, jj) * p(gt,w,h, ii, jj);
	}
}

static void compute_structure_tensor_field(float *st,
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, int w, int h)
{

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float atwa[STLEN];
		compute_structure_tensor_here(atwa, wv, wo, kside, gx, gy, w, h,
				i, j);
		for (int l = 0; l < STLEN; l++)
			setpixel_float_image_vec(st, w, h, STLEN, i, j, l, atwa[l]);
	}
}

static void compute_rhs(float *rhs,
		float *wv, int (*wo)[2], int kside,
		float *gx, float *gy, float *gt, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float rhsh[2];
		compute_rhs_here(rhsh, wv, wo, kside, gx, gy, gt, w, h,
				i, j);
		setpixel_float_image_vec(rhs, w, h, 2, i, j, 0, rhsh[0]);
		setpixel_float_image_vec(rhs, w, h, 2, i, j, 1, rhsh[1]);
	}
}

static void solve_pointwise(float *u, float *v, float *st, float *rhs,
		int w, int h)
{
	float (*x_st)[w][STLEN] = (void*)st;
	float (*x_rhs)[w][2] = (void*)rhs;
	float (*x_u)[w] = (void*)u;
	float (*x_v)[w] = (void*)v;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float f[2];
		solve_sdp_2x2(f, x_st[j][i], x_rhs[j][i]);
		x_u[j][i] = f[0];//x_st[j][i][0];
		x_v[j][i] = f[1];//x_st[j][i][1];
	}
}

static void least_squares_ofc(float *u, float *v,
		float *a, float *b, int w, int h,
		int kside, float sigma)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	if (kside == -1) { // global approximation

	} else {
		if (kside % 2 != 1) exit(fprintf(stderr,
			"I need an ODD window size (got %d)\n", kside));
		if (kside > 40) exit(fprintf(stderr,
			"window size %d too large\n", kside));
		int wo[kside*kside][2]; fill_window_offsets(wo, kside);
		float wv[kside*kside]; fill_window_values(wv, wo, kside, sigma);
		float *st = xmalloc(w * h * STLEN * sizeof(float));
		float *rhs = xmalloc(w * h * 2 * sizeof(float));
		compute_structure_tensor_field(st, wv, wo, kside, gx, gy, w, h);
		//iio_save_image_float_vec("/tmp/st.tiff", st, w, h, STLEN);
		compute_rhs(rhs, wv, wo, kside, gx, gy, gt, w, h);
		//iio_save_image_float_vec("/tmp/rhs.tiff", rhs, w, h, 2);
		solve_pointwise(u, v, st, rhs, w, h);
		//iio_save_image_float("/tmp/u.tiff", u, w, h);
		//iio_save_image_float("/tmp/v.tiff", v, w, h);
		free(rhs);
		free(st);
	}
	free(gx);
	free(gy);
       	free(gt);
}

int main(int argc, char *argv[])
{
	if (argc != 6)
		exit(fprintf(stderr, "usage:\n\t%s kside sigma a b f\n",*argv));
	int kside = atoi(argv[1]);
	float sigma = atoi(argv[2]);
	char *filename_a = argv[3];
	char *filename_b = argv[4];
	char *filename_f = argv[5];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		exit(fprintf(stderr, "input images size mismatch\n"));
	iio_save_image_float("/tmp/merdota.tiff", a, w, h);
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	least_squares_ofc(u, v, a, b, w, h, kside, sigma);
	float *f = xmalloc(w * h * 2 * sizeof(float));
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_save_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
