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
static void compute_input_derivatives(float *gx, float *gy, float *gt,
		float *a, float *b, int w, int h)
{
	float (*Ex)[w] = (void*)gx;
	float (*Ey)[w] = (void*)gy;
	float (*Et)[w] = (void*)gt;
	extension_operator_float p = extend_float_image_by_zero;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		Ey[j][i] = (1.0/4) * (
				p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j)
				);
		Ex[j][i] = (1.0/4) * (
				p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1)
				);
		Et[j][i] = (1.0/4) * (
				p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1)
				);
	}
}
static void compute_bar(float *ubar, float *u, int w, int h)
{
	float (*ub)[w] = (void*)ubar;
	extension_operator_float p = extend_float_image_constant;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ub[j][i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}
static void hs_iteration(float *u, float *v, float *a, float *b,
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
	free(ubar); free(vbar);
}
static void hs(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (int i = 0; i < w*h; i++) u[i] = v[i] = 0;
	for (int i = 0; i < niter; i++) hs_iteration(u, v, a, b, gx, gy, gt, w, h, alpha);
	free(gx); free(gy); free(gt);
}
