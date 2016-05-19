
float rip_product(float *g, int w, int h, float x[2], float p[2], float q[2])
{
	// observation: piecewise constant interpolation
	int i = lrint(x[0]);
	int j = lrint(x[1]);
	float E = g[3*(j*w + i) + 0];
	float F = g[3*(j*w + i) + 1];
	float G = g[3*(j*w + i) + 2];
	return E * p[0]*q[0] + F * (p[1]*q[0] + p[0]*q[1]) + G * p[1]*q[1];
}

float rip_curve_length(float *g, int w, int h, float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n-1; i++)
	{
		float *xi = x + 2*i;
		float *xI = x + 2*(i+1);
		float vi[2] = { xI[0] - xi[0], xI[1] - xi[1] };
		float dt = hypot(vi[0], vi[1]);
		r += sqrt(rip_product(g, w, h, xi, vi, vi)) / dt;
	}
	return r;
}

float rip_curve_energy(float *g, int w, int h, float *x, int n)
{
	long double r = 0;
	for (int i = 0; i < n-1; i++)
	{
		float *xi = x + 2*i;
		float *xI = x + 2*(i+1);
		float vi[2] = { xI[0] - xi[0], xI[1] - xi[1] };
		float dt = hypot(vi[0], vi[1]);
		r += rip_product(g, w, h, xi, vi, vi) / dt;
	}
	return r;
}
