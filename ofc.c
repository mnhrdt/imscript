// implementation and tests of the optical flow constraint

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"


#include "fragments.c"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

// temporal derivative, forward differences
static void fill_it(float *it, float *xa, float *xb, int w, int h)
{
	FORI(w*h)
		it[i] = xb[i] - xa[i];
}


typedef float (*extension_operator_float)(float*,int,int,int,int);

static float extend_float_image_by_zero(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0 || j < 0 || i > w-1 || j > h-1)
		return 0;
	else
		return x[j][i];
}

static float extend_float_image_constant(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j][i];
}


// spatial derivatives, centered differences
static void fill_gradx(float *out_gradx, float *in_x, int w, int h)
{
	float (*gradx)[w][2] = (void*)out_gradx;
	//float (*x)[w] = (void*)in_x;
	float *x = in_x;

	extension_operator_float p = extend_float_image_by_zero;

	FORJ(h) FORI(w) {
		gradx[j][i][0] = (p(x,w,h, i+1, j) - p(x,w,h, i-1, j))/2;
		gradx[j][i][1] = (p(x,w,h, i, j+1) - p(x,w,h, i, j-1))/2;
	}
}

static void fill_git(float *out_git, float *xa, float *xb, int w, int h)
{
	float (*git)[w][2] = (void*)out_git;
	float (*grad0)[w][2] = xmalloc(w*h*2*sizeof(float));
	float (*grad1)[w][2] = xmalloc(w*h*2*sizeof(float));
	fill_gradx(grad0[0][0], xa, w, h);
	fill_gradx(grad1[0][0], xb, w, h);
	FORJ(h) FORI(w) FORL(2)
		git[j][i][l] = grad1[j][i][l] - grad0[j][i][l];
	xfree(grad0);
	xfree(grad1);
}


static void fill_hess(void *out_hess, void *in_x, int w, int h)
{
	float (*hess)[w][3] = out_hess;
	//float (*x)[w][2] = in_x;
	float *x = in_x;

	extension_operator_float p = extend_float_image_by_zero;

	FORJ(h) FORI(w) {
		hess[j][i][0] = p(x,w,h, i-1, j)
                           -2 * p(x,w,h, i, j)
                              + p(x,w,h, i+1, j);
		hess[j][i][3] = p(x,w,h, i, j-1)
                           -2 * p(x,w,h, i, j)
                              + p(x,w,h, i, j+1);
		hess[j][i][2] = 0.25 * (
				p(x,w,h, i-1, j-1) + p(x,w,h, i+1, j+1) -
				p(x,w,h, i+1, j-1) - p(x,w,h, i-1, j+1) ); 
	}
}

#include "vvector.h"

static float solve_linear_2x2(float x[2], float a[2][2], float b[2])
{
	float inva[2][2], det;
	INVERT_2X2(inva, det, a);
	x[0] = inva[0][0] * b[0] + inva[0][1] * b[1];
	x[1] = inva[1][0] * b[0] + inva[1][1] * b[1];
	return det;
}

static void fill_fbar(float *out_ubar, float *u, int w, int h)
{
	float (*ubar)[w] = (void*)out_ubar;

	extension_operator_float p = extend_float_image_constant;

	FORJ(h) FORI(w)
		ubar[j][i] = (1.0/6) * (
				p(u,w,h, i-1, j)+
				p(u,w,h, i, j-1)+
				p(u,w,h, i+1, j)+
				p(u,w,h, i, j+1)
			) + (1.0/12) * (
				p(u,w,h, i-1, j-1)+
				p(u,w,h, i+1, j-1)+
				p(u,w,h, i+1, j+1)+
				p(u,w,h, i-1, j+1)
			);
}

static void hs_iteration(float *out_nextu, float *out_nextv,
		float *in_u, float *in_v,
		float *in_ex, float *in_ey, float *in_et,
		int w, int h, float alpha)
{
	float *ubar = xmalloc(w*h*sizeof*ubar);
	float *vbar = xmalloc(w*h*sizeof*vbar);
	fill_fbar(ubar, in_u, w, h);
	fill_fbar(vbar, in_v, w, h);

	FORJ(h) FORI(w) {
		int idx = j*w + i;
		float ex = in_ex[idx];
		float ey = in_ey[idx];
		float et = in_et[idx];
		float ub = ubar[idx];
		float vb = vbar[idx];
		float fac = alpha*alpha + ex*ex + ey*ey;
		float nextu = ub - ex * (ex * ub + ey * vb + et)/fac;
		float nextv = vb - ey * (ex * ub + ey * vb + et)/fac;
		out_nextu[idx] = nextu;
		out_nextv[idx] = nextv;
	}

	free(ubar);
	free(vbar);
}

static void ofc(float *xa, float *xb, int w, int h)
{
	float (*x)[w] = (void*)xa;
	float (*y)[w] = (void*)xb;

	float (*it)[w] = xmalloc(w*h*sizeof(float));
	float (*gradx)[w][2] = xmalloc(w*h*2*sizeof(float));
	float (*hessx)[w][3] = xmalloc(w*h*3*sizeof(float);
	float (*gitx)[w][2] = xmalloc(w*h*2*sizeof(float));
	float (*f1)[w][2] = xmalloc(w*h*2*sizeof(float));
	float (*f2)[w][2] = xmalloc(w*h*2*sizeof(float));

	float (*u)[w] = xmalloc(w*h*sizeof(float));
	float (*v)[w] = xmalloc(w*h*sizeof(float));
	float (*ubar)[w] = xmalloc(w*h*sizeof(float));
	float (*vbar)[w] = xmalloc(w*h*sizeof(float));

	fill_it(it[0], xa, xb, w, h);
	fill_gradx(gradx[0][0], xa, w, h);
	fill_hess(hessx[0][0], xa, w, h);
	fill_git(gitx[0][0], xa, xb, w, h);

	// compute f1
	FORJ(h) FORI(w) {
		float k = -it[j][i];
		float a = gradx[j][i][0];
		float b = gradx[j][i][1];
		f1[j][i][0] = (-k * a)/(a*a + b*b);
		f1[j][i][1] = (-k * b)/(a*a + b*b);
	}

	// compute f2
	FORJ(h) FORI(w) {
		float a[2][2] = {{hessx[j][i][0], hessx[j][i][1]},
			{hessx[j][i][1], hessx[j][i][2]}};
		float b[2] = {-gitx[j][i][0], -gitx[j][i][1]};
		float res[2], det;
		det = solve_linear_2x2(res, a, b);
		FORL(2)
			f2[j][i][l] = res[l];
	}

	iio_save_image_float_vec("/tmp/it.tiff", it[0], w, h, 1);
	iio_save_image_float_vec("/tmp/grad.tiff", gradx[0][0], w, h, 2);
	iio_save_image_float_vec("/tmp/f1.tiff", f1[0][0], w, h, 2);
	iio_save_image_float_vec("/tmp/f2.tiff", f2[0][0], w, h, 2);

	xfree(it);
	xfree(gradx);
	xfree(hessx);
	xfree(gitx);
	xfree(f1);
	xfree(f2);
}

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s framea frameb\n", *v);
		return EXIT_FAILURE;
	}
	char *framea = v[1];
	char *frameb = v[2];

	int w[2], h[2], pd[2];
	float *xa = iio_read_image_float_vec(framea, w, h, pd);
	float *xb = iio_read_image_float_vec(frameb, w+1, h+1, pd+1);

	if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
		error("two frames sizing mismatch");

	if (*pd != 1) error("only gray level by now");

	ofc(xa, xb, *w, *h);

	return EXIT_SUCCESS;
}
