#include <assert.h>
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

	if (0) {
		FILE *f = fopen("/tmp/lk.kkk", "w");
		fprintf(f, "P2\n%d %d\n65535\n", kside, kside);
		for (int i = 0; i < n; i++)
			fprintf(f, "%g\n", wv[i]);
		fclose(f);
	}
}

#include "svd.c"
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


/*  solvps.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
int solvps(double *a,double *b,int n)
{ double *p,*q,*r,*s,t;
  int j,k;
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(q=a+j*n; q<p ;++q) *p-= *q* *q;
    if(*p<=0.) return -1;
    *p=sqrt(*p);
    for(k=j+1,q=p+n; k<n ;++k,q+=n){
      for(r=a+j*n,s=a+k*n,t=0.; r<p ;) t+= *r++ * *s++;
      *q-=t; *q/= *p;
     }
   }
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(k=0,q=a+j*n; k<j ;) b[j]-=b[k++]* *q++;
    b[j]/= *p;
   }
  for(j=n-1,p=a+n*n-1; j>=0 ;--j,p-=n+1){
    for(k=j+1,q=p+n; k<n ;q+=n) b[j]-=b[k++]* *q;
    b[j]/= *p;
   }
  return 0;
}

static double solve_sdp_6x6(double x[6], double A[6][6], double b[6])
{
	fprintf(stderr,"A = \n"); for (int j = 0; j < 6; j++) {
		for (int i = 0; i < 6; i++)fprintf(stderr," %g", A[j][i]);fprintf(stderr,"\n");}
	fprintf(stderr,"rhs = \n");for(int i=0;i<6;i++)fprintf(stderr," %g",b[i]);fprintf(stderr,"\n");

	int r = solvps(A[0], b, 6);
	if(r<0)exit(fprintf(stderr,"affine structure tensor is singular\n"));
	for (int i = 0; i < 6; i++)
		x[i] = b[i];
	//double d[6], u[6][6], v[6][6];
	//svd(d, A[0], u[0], 6, v[0], 6);
	//printf("u = \n"); for (int j = 0; j < 6; j++) {
	//	for (int i = 0; i < 6; i++)printf(" %g", u[j][i]);printf("\n");}
	//printf("v = \n"); for (int j = 0; j < 6; j++) {
	//	for (int i = 0; i < 6; i++)printf(" %g", v[j][i]);printf("\n");}
	//printf("d = \n");for(int i=0;i<6;i++)printf(" %g",d[i]);printf("\n");
	fprintf(stderr,"x = \n");for(int i=0;i<6;i++)fprintf(stderr," %g",x[i]);fprintf(stderr,"\n");
	return 0;
}

static double solve_sdp_nxn(double *x, double *A, double *b, int n)
{
	//fprintf(stderr, "\n\nSOLVE SDP NXN (n=%d)\n\n", n);
	//fprintf(stderr,"A = \n"); for (int j = 0; j < n; j++) {
	//	for (int i = 0; i < n; i++)fprintf(stderr," %g", A[j*n+i]);fprintf(stderr,"\n");}
	//fprintf(stderr,"rhs = \n");for(int i=0;i<n;i++)fprintf(stderr," %g",b[i]);fprintf(stderr,"\n");

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			assert(A[n*i+j] == A[n*j + i]);

	int r = solvps(A, b, n);
	for (int i = 0; i < n; i++)
		x[i] = b[i];

	double d[n], u[n][n], v[n][n];
	svd(d, A, u[0], n, v[0], n);
	//fprintf(stderr,"u = \n"); for (int j = 0; j < n; j++) {
	//	for (int i = 0; i < n; i++)fprintf(stderr," %g", u[j][i]);fprintf(stderr,"\n");}
	//fprintf(stderr,"v = \n"); for (int j = 0; j < n; j++) {
	//	for (int i = 0; i < n; i++)fprintf(stderr," %g", v[j][i]);fprintf(stderr,"\n");}
	//fprintf(stderr,"d = \n");for(int i=0;i<n;i++)fprintf(stderr," %g",d[i]);fprintf(stderr,"\n");
	//fprintf(stderr,"x = \n");for(int i=0;i<n;i++)fprintf(stderr," %g",x[i]);fprintf(stderr,"\n");
	if(r<0)exit(fprintf(stderr,"polynomial structure tensor is singular\n"));
	return 0;
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
#pragma omp parallel for
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
#pragma omp parallel for
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

#pragma omp parallel
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float f[2];
		solve_sdp_2x2(f, x_st[j][i], x_rhs[j][i]);
		x_u[j][i] = f[0];//x_st[j][i][0];
		x_v[j][i] = f[1];//x_st[j][i][1];
	}
}

static void global_constant_approximation(float *u, float *v,
		float *gx, float *gy, float *gt, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;

	float atwa[STLEN], rhs[2];
	atwa[0] = atwa[1] = atwa[2] = rhs[0] = rhs[1] = 0;

	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			atwa[0] += 1 * sqr(p(gx, w, h, i, j));
			atwa[1] += 1 * p(gx,w,h, i, j) * p(gy,w,h, i, j);
			atwa[2] += 1 * sqr(p(gy, w, h, i, j));
			rhs[0] -= 1 * p(gx,w,h, i, j) * p(gt,w,h, i, j);
			rhs[1] -= 1 * p(gy,w,h, i, j) * p(gt,w,h, i, j);
		}

	float f[2];
	solve_sdp_2x2(f, atwa, rhs);
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			u[j*w + i] = f[0];
			v[j*w + i] = f[1];
		}
}

static void global_affine_approximation(float *u, float *v,
		float *gx, float *gy, float *gt, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;

	double ast[6][6], rhs[6] = {0,0,0,0,0,0};
	for (int j = 0; j < 6; j++)
	for (int i = 0; i < 6; i++)
		ast[j][i] = 0;

	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			double x = i;
			double y = j;
			double Ex = p(gx, w, h, i, j);
			double Ey = p(gy, w, h, i, j);
			double Et = p(gt, w, h, i, j);
			ast[0][0] += Ex * Ex * x * x;
			ast[0][1] += Ex * Ex * x * y;
			ast[0][2] += Ex * Ex * x;
			ast[0][3] += Ex * Ey * x * x;
			ast[0][4] += Ex * Ey * x * y;
			ast[0][5] += Ex * Ey * x;
			ast[1][0] += Ex * Ex * x * y;
			ast[1][1] += Ex * Ex * y * y;
			ast[1][2] += Ex * Ex * y;
			ast[1][3] += Ex * Ey * x * y;
			ast[1][4] += Ex * Ey * y * y;
			ast[1][5] += Ex * Ey * y;
			ast[2][0] += Ex * Ex * x;
			ast[2][1] += Ex * Ex * y;
			ast[2][2] += Ex * Ex;
			ast[2][3] += Ex * Ey * x;
			ast[2][4] += Ex * Ey * y;
			ast[2][5] += Ex * Ey;
			ast[3][0] += Ex * Ey * x * x;
			ast[3][1] += Ex * Ey * x * y;
			ast[3][2] += Ex * Ey * x;
			ast[3][3] += Ey * Ey * x * x;
			ast[3][4] += Ey * Ey * x * y;
			ast[3][5] += Ey * Ey * x;
			ast[4][0] += Ex * Ey * x * y;
			ast[4][1] += Ex * Ey * y * y;
			ast[4][2] += Ex * Ey * y;
			ast[4][3] += Ey * Ey * x * y;
			ast[4][4] += Ey * Ey * y * y;
			ast[4][5] += Ey * Ey * y;
			ast[5][0] += Ex * Ey * x;
			ast[5][1] += Ex * Ey * y;
			ast[5][2] += Ex * Ey;
			ast[5][3] += Ey * Ey * x;
			ast[5][4] += Ey * Ey * y;
			ast[5][5] += Ey * Ey;
			rhs[0] -= Et * Ex * x;
			rhs[1] -= Et * Ex * y;
			rhs[2] -= Et * Ex;
			rhs[3] -= Et * Ey * x;
			rhs[4] -= Et * Ey * y;
			rhs[5] -= Et * Ey;
		}

	//for (int j = 0; j < 6; j++)
	//for (int i = 0; i < 6; i++)
	//	ast[j][i] /= sqrt(w*h);
	//for (int i = 0; i < 6; i++)
	//	rhs[i] /= sqrt(w*h);

	double f[6];
	solve_sdp_6x6(f, ast, rhs);
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			double x = i;
			double y = j;
			u[j*w + i] = f[0]*x + f[1]*y + f[2];
			v[j*w + i] = f[3]*x + f[4]*y + f[5];
		}
}


static double ipow(double base, int exponent)
{
	double result = 1;
	while (exponent)
	{
		if (exponent & 1)
			result *= base;
		exponent >>= 1;
		base *= base;
	}
	return result;
}

static double monomium(double x, double y, int exponent[2])
{
	return ipow(x,exponent[0]) * ipow(y,exponent[1]);
}

static double polynomium(double x, double y,
		double *f, int (*polindex)[2], int nc)
{
	double r = 0;
	for (int i = 0; i < nc; i++)
		r += f[i] * monomium(x, y, polindex[i]);
	return r;
}

static void global_polynomial_approximation(float *u, float *v,
		float *gx, float *gy, float *gt, int w, int h, int deg)
{
	int nc = (deg + 2) * (deg + 1) / 2; // number of coefficients
	double pst[2*nc][2*nc], rhs[2*nc];
	int polindex[nc][2], idx = 0;
	for (int j = 0; j <= deg; j++)
	for (int i = 0; i <= j; i++)
	{
		polindex[idx][0] = j-i;
		polindex[idx][1] = i;
		idx += 1;
	}
	assert(idx == nc);

	//for (int i = 0; i < nc; i++)
	//	fprintf(stderr, "polindex[%d] = {%d, %d};\n",
	//			i, polindex[i][0], polindex[i][1]);

	for (int j = 0; j < 2*nc; j++)
	for (int i = 0; i < 2*nc; i++)
		pst[j][i] = 0;
	for (int i = 0; i < 2*nc; i++)
		rhs[i] = 0;

	double nfac = sqrt(w*h);//100;

	int passepartoutw = 2;//0.004 * w;
	int passepartouth = 2;//0.004 * h;

	for (int jj = passepartouth; jj < h - passepartouth; jj++)
		for (int ii = passepartoutw; ii < w - passepartoutw; ii++)
		{
			extension_operator_float p= extend_float_image_constant;
			double x = ii/nfac;
			double y = jj/nfac;
			double Ex = p(gx, w, h, ii, jj);
			double Ey = p(gy, w, h, ii, jj);
			double Et = p(gt, w, h, ii, jj);

			for (int j = 0; j < nc; j++)
			for (int i = 0; i < nc; i++)
			{
				int m[2] = {
					polindex[i][0] + polindex[j][0],
					polindex[i][1] + polindex[j][1]
				};
				pst[j][i] += Ex * Ex * monomium(x, y, m);
				pst[j+nc][i] += Ex * Ey * monomium(x, y, m);
				pst[j][i+nc] += Ex * Ey * monomium(x, y, m);
				pst[j+nc][i+nc] += Ey * Ey * monomium(x, y, m);
			}
			for (int i = 0; i < nc; i++)
			{
				int m[2] = { polindex[i][0], polindex[i][1] };
				rhs[i] -= Et * Ex * monomium(x, y, m);
				rhs[i+nc] -= Et * Ey * monomium(x, y, m);
			}
		}

	double f[2*nc];
	solve_sdp_nxn(f, pst[0], rhs, 2*nc);
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			double x = i/nfac;
			double y = j/nfac;
			u[j*w + i] = polynomium(x, y, f, polindex, nc);;
			v[j*w + i] = polynomium(x, y, f+nc, polindex, nc);;
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
	if (kside == -1)
		global_constant_approximation(u, v, gx, gy, gt, w, h);
	else if (kside == -2)
		global_affine_approximation(u, v, gx, gy, gt, w, h);
	else if (kside == -3) {
		int deg = sigma;
		fprintf(stderr, "USING POLYNOMIALS OF DEGREE %d\n", deg);
		global_polynomial_approximation(u, v, gx, gy, gt, w, h, deg);
	}
	else {
		if (kside % 2 != 1) exit(fprintf(stderr,
			"I need an ODD window size (got %d)\n", kside));
		if (kside > 70) exit(fprintf(stderr,
			"window size %d too large\n", kside));
		int wo[kside*kside][2]; fill_window_offsets(wo, kside);
		float wv[kside*kside]; fill_window_values(wv, wo, kside, sigma);
		float *st = xmalloc(w * h * STLEN * sizeof(float));
		float *rhs = xmalloc(w * h * 2 * sizeof(float));
		compute_structure_tensor_field(st, wv, wo, kside, gx, gy, w, h);
		iio_write_image_float_vec("stf.tiff", st, w, h, STLEN);
		compute_rhs(rhs, wv, wo, kside, gx, gy, gt, w, h);
		iio_write_image_float_vec("rhs.tiff", rhs, w, h, 2);
		solve_pointwise(u, v, st, rhs, w, h);
		//iio_write_image_float("/tmp/u.tiff", u, w, h);
		//iio_write_image_float("/tmp/v.tiff", v, w, h);
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
	//iio_write_image_float("/tmp/merdota.tiff", a, w, h);
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));
	least_squares_ofc(u, v, a, b, w, h, kside, sigma);
	float *f = xmalloc(w * h * 2 * sizeof(float));
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_write_image_float_vec(filename_f, f, w, h, 2);
	return EXIT_SUCCESS;
}
