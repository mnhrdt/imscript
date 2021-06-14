// implementation of Cohen-Lee-Miller-Pachocki-Sidford geometric medians
// ref. "Geometric Median in Nearly Linear Time", arXiv:160605225v

#include <assert.h>
#include <math.h>

#include "linalg.c" // cholesky, solve_spd


static void linear_average(float *o, int d, int n, float a[n][d])
{
	for (int k = 0; k < d; k++)
	{
		o[k] = 0;
		for (int i = 0; i < n; i++)
			o[k] += a[i][k] / n;
	}
}

// y[k] = (1/n) * sum_i x[i][k]
static void float_avg(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

static float fnorm(float *x, int n)
{
	switch(n) {
	case 1: return fabs(x[0]);
	case 2: return hypot(x[0], x[1]);
	default: return hypot(x[0], fnorm(x+1, n-1));
	}
}

static float medscore(float *xx, int idx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	float r = 0;
	for (int i = 0; i < n; i++)
		if (i != idx)
		{
			float v[d];
			for (int j = 0; j < d; j++)
				v[j] = x[idx][j] - x[i][j];
			r += fnorm(v, d);
		}
	return r;
}

// y[] = x[i][] which is closest to the euclidean median
// runs in quadratic time
static void float_med(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	int midx = 0;
	float misc = medscore(xx, 0, d, n);
	for (int i = 1; i < n; i++)
	{
		float si = medscore(xx, i, d, n);
		if (si < misc) {
			midx = i;
			misc = si;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];

}

static float euclidean_distance(float *p, float *q, int d)
{
	double r = 0;
	for (int k = 0; k < d; k++)
		r = hypot(r, p[k] - q[k]);
	return r;
}

static float objective_function(int d, int n, float a[n][d], float *x)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += euclidean_distance(a[i], x, d);
	return r;
}

static float penalized_objective(float *a[], int d, int n, float *x,
		float e, float t)
{
	double r = 0;
	for (int i = 0; i < n; i++)
	{
		float e = euclidean_distance(x, a[i], d);
		float s = sqrt(1 + t*t * e*e);
		r += s - log(1 + s);
	}
	return r;
}

void clmps_median(
		float *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		float a[n][d],    // input list of points
		float e           // desired accuracy 0<e<1
		)
{
	assert(0 < e && e < 1);

	// compute a 2-approximate geometric median and use it to center
	float x[d];
	linear_average(x, d, n, a);
	float f_star = objective_function(d, n, a, x);
	float t[n];
	for (int i = 0; i < n; i++)
		t[i] = pow(601.0/600, i)/(400*f_star);
	float e_star = e/3;
	float t_star = 2 * n / (e_star * f_star);
	float e_v = pow(e_star/(7*n), 2)/8;
	float e_c = pow(e_v/36, 1.5);
	float y[d];
	//line_search(x, t[0], t[0], 0, e_c);


}

#include "random.c"
#include <stdio.h>

static void get_random_subset(int *s, int n, int k)
{
	// knuth O(n) naive method
	for (int i = 0; i < n; i++)
		if (randombounds(0, n-i-1) < k)
			s[--k] = i;
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

float crude_approximate(
		float *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		float a[n][d],    // input list of points
		int K             // subset size
		)
{
	assert(K > 0 && K < n);

	int S1[K], S2[K];
	get_random_subset(S1, n, K);
	get_random_subset(S2, n, K);

	float α[K][K];
	for (int i = 0; i < K; i++)
	for (int j = 0; j < K; j++)
		α[i][j] = euclidean_distance(a[S2[i]], a[S1[j]], d);

	for (int i = 0; i < K; i++)
		qsort(α[i], K, sizeof*α[i], compare_floats);

	int p = lrint(0.65*K);
	int ι = 0;
	for (int i = 0; i < K; i++)
		if (α[i][p] < α[ι][p])
			ι = i;

	for (int i = 0; i < d; i++)
		out[i] = a[S2[ι]][i];
	return α[ι][p];
}

void linear_approximate_median(
		float *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		float a[n][d],    // input list of points
		float ε           // desired precision
		)
{
	float T = pow(60/ε, 2);
	float x[d];
	float λ = crude_approximate(x, d, n, a, lrint(sqrt(T)));
	float η = 6*λ*sqrt(2/T)/n;
	{//debug
		fprintf(stderr, "\td = %d\n", d);
		fprintf(stderr, "\tn = %d\n", n);
		fprintf(stderr, "\tε = %g\n", ε);
		fprintf(stderr, "\tT = %g\n", T);
		fprintf(stderr, "\trT = %ld\n", lrint(sqrt(T)));
		fprintf(stderr, "linear approximate median\n");
		fprintf(stderr, "\tλ = %g\n", λ);
		fprintf(stderr, "\tη = %g\n", η);
		fprintf(stderr, "\tR = %g\n", 6*λ);
		if (d == 2)
			fprintf(stderr, "\tx0 = %g %g\n", x[0], x[1]);
	}

	for (int k = 0; k < T; k++)
	{
		int i = randombounds(0, n-1);
		float ρ = euclidean_distance(x, a[i], d);
		float g[d];
		for (int j = 0; j < d; j++)
			g[j] = ρ>0 ? n * (x[j] - a[i][j])/ρ : 0;
		for (int j = 0; j < d; j++)
			x[j] = x[j] - η*g[j];
		// TODO: add reprojection to the sphere of radious 6λ
		if (d==2 && (k<10 || k+10>T))
			fprintf(stderr, "\tx_%d = %lf %lf\n", k, x[0], x[1]);
	}

	for (int j = 0; j < d; j++)
		out[j] = x[j];
}

// euclidean distance between the vectors x and y, regularized around 0
static float fdiste(float *x, float *y, int n, float e)
{
	return n ? hypot(*x - *y, fdiste(x + 1, y + 1, n - 1, e)) : e;
}

#include "smapa.h"
SMART_PARAMETER_SILENT(WEISZ_NITER,10)

// y[k] = euclidean median of the vectors x[i][k]
static void float_weisz(float *y, float *x, int d, int n, int N, float ε)
{
	float_avg(y, x, d, n);
	for (int k = 0; k < N; k++) {
		float a[d], b = 0;
		for (int l = 0; l < d; l++)
			a[l] = 0;
		for (int i = 0; i < n; i++) {
			float dxy = fdiste(x + i*d, y, d, ε);
			for (int l = 0; l < d; l++)
				a[l] += x[i*d + l]/dxy;
			b += 1/dxy;
		}
		for (int l = 0; l < d; l++)
			y[l] = a[l]/b;
		if (d == 2 && N < 100)
			fprintf(stderr, "w[%d] = %lf %lf\n", k, y[0], y[1]);
	}
}

static void fill_points(int d, int n, float x[n][d])
{
	assert(d == 2);

	int m = 3;         // number of centers
	float c[][2] = {  // centers
		{100, 240},
		{400, 250},
		{100, 260},
	};
	float σ = 4;

	for (int i = 0; i < n; i++)
	{
		int j = i % m;
		assert(j >= 0 && j < m);
		x[i][0] = c[j][0] + σ * random_normal();
		x[i][1] = c[j][1] + σ * random_normal();
	}

	FILE *f = fopen("/tmp/puns.txt", "w");
	for (int i = 0; i < n; i++)
		fprintf(f, "%lf %lf\n", x[i][0], x[i][1]);
	fclose(f);
}

static void fill_points2(int d, int n, float x[n][d])
{
	assert(d == 2);

	int m = 2;         // number of centers
	float c[][2] = {  // centers
		{100, 100},
		{400, 400}
	};
	float σ[] = {15, 3};

	for (int i = 0; i < n; i++)
	{
		int j = i % m;
		assert(j >= 0 && j < m);
		x[i][0] = c[j][0] + σ[j] * random_normal();
		x[i][1] = c[j][1] + σ[j] * random_normal();
	}

	FILE *f = fopen("/tmp/puns2.txt", "w");
	for (int i = 0; i < n; i++)
		fprintf(f, "%lf %lf\n", x[i][0], x[i][1]);
	fclose(f);
}

#include "pickopt.c"
int main_clmps(int c, char *v[])
{
	float weisz_n = atoi(pick_option(&c, &v, "wn", "10"));
	float weisz_e = atof(pick_option(&c, &v, "we", "0.01"));
	float crude_K = atoi(pick_option(&c, &v, "ck", "10"));
	float linap_e = atof(pick_option(&c, &v, "le", "0.5"));

	int d = 2;     // dimension
	int n = 3*10000; // total number of points
	float x[n][d];
	fill_points2(d, n, x);

	float E;

	float avg[d];
	linear_average(avg, d, n, x);
	E = objective_function(d, n, x, avg)/n;
	fprintf(stderr, "avg = %lf %lf (%lf)\n", avg[0], avg[1], E);

	float wei[d];
	float_weisz(wei, *x, d, n, weisz_n, weisz_e);
	E = objective_function(d, n, x, wei)/n;
	fprintf(stderr, "wei = %lf %lf (%lf)\n", wei[0], wei[1], E);

	float cru[d];
	float r = crude_approximate(cru, d, n, x, crude_K);
	E = objective_function(d, n, x, cru)/n;
	fprintf(stderr, "cru = %lf %lf {%lf} (%lf)\n", cru[0], cru[1], r, E);

	float lin[d];
	linear_approximate_median(lin, d, n, x, linap_e);
	E = objective_function(d, n, x, lin)/n;
	fprintf(stderr, "lin = %lf %lf (%lf)\n", lin[0], lin[1], E);

	float med[d];
	float_med(med, *x, d, n);
	E = objective_function(d, n, x, med)/n;
	fprintf(stderr, "med = %lf %lf (%lf)\n", med[0], med[1], E);


	return 0;
}

//static float objective_function(int d, int n, float a[n][d], float *x)
//{
//	double r = 0;
//	for (int i = 0; i < n; i++)
//		r += euclidean_distance(a[i], x, d);
//	return r;
//}

static void gradient(float *g, int d, int n, float a[n][d], float *x)
{
	for (int k = 0; k < d; k++)
		g[k] = 0;
	for (int i = 0; i < n; i++)
	{
		double r = 0;
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);
		for (int k = 0; k < d; k++)
			g[k] += (x[k] - a[i][k]) / r;
	}
}

static float weiszfeld_weight(int d, int n, float a[n][d], float *x)
{
	double R = 0; // sum of all inverse norms x-ai
	for (int i = 0; i < n; i++)
	{
		double r = 0; // norm of x-ai
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);
		R += 1/r;
	}
	return R;
}

static void hessian(float *h, int d, int n, float a[n][d], float *x)
{
	float (*H)[d] = (void*)h;

	// initialize to the identity
	for (int k = 0; k < d; k++)
	for (int l = 0; l < d; l++)
		H[k][l] = k == l;

	// first term: a scalar multiple of the identity
	double R = 0; // sum of all inverse norms x-ai
	for (int i = 0; i < n; i++)
	{
		double r = 0; // norm of x-ai
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);
		R += 1/r;
	}
	for (int k = 0; k < d; k++)
		H[k][k] *= R;

	if (d == 2)
		fprintf(stderr, "\t\tH=%g %g %g %g\n", h[0],h[1],h[2],h[3]);

	// second-term: accumulate singular perturbations
	for (int i = 0; i < n; i++)
	{
		double r = 0; // norm of x-ai (TODO: do not recompute)
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);

		for (int k = 0; k < d; k++)
		for (int l = 0; l < d; l++)
			H[k][l] -= (x[k]-a[i][k])*(x[l]-a[i][l]) / (r*r*r);
	}

	if (d == 2)
		fprintf(stderr, "\t\tH=%g %g %g %g\n", h[0],h[1],h[2],h[3]);
}

#include "iio.h"

#include "smapa.h"
SMART_PARAMETER_SILENT(XMIN,-2)
SMART_PARAMETER_SILENT(XMAX,12)
SMART_PARAMETER_SILENT(YMIN,-2)
SMART_PARAMETER_SILENT(YMAX,12)
SMART_PARAMETER_SILENT(NUMIT,20)

// fancier weiszfeld variands based on gradient descent
// (by default, plain weiszfeld)
//
// Variants:
// 	stepsize-constant: use a constant stepsize (relative to the gradient)
// 	stepsize-absolute: use a constant-length stepsize (non-convergent!)
// 	stepsize-factor: multiply the Weiszfeld lambda by this factor
// 	search-armijo: linear armijo search with parameters = 1/2
// 	stochastic-computations: 
//
int main_weisz(int c, char *v[])
{
	char *out_sampling = pick_option(&c, &v, "o", "");

	if (c != 1)
		return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);

	int n, d;
	void *aa = iio_read_image_float("-", &d, &n);
	float (*a)[d] = aa;

	if (d == 2 && *out_sampling)
	{
		int w = 1000;
		int h = 1000;
		float *o = malloc(w*h*sizeof*o);
		for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			float x[2] = {
				XMIN() + i*(XMAX()-XMIN())/w,
				YMIN() + j*(YMAX()-YMIN())/h,
			};
			o[j*w+i] = objective_function(d, n, a, x);
		}
		iio_write_image_float(out_sampling, o, w, h);
		free(o);
	}

	float x[d]; // initialization
	for (int i = 0; i < d; i++)
		x[i] = 0;

	float E = objective_function(d, n, a, x);

	fprintf(stderr, "got %d points in dimension %d\n", n, d);
	fprintf(stderr, "energy at zero = %g\n", E);


	if (n < 10)
		for (int i = 0; i < n; i++)
			fprintf(stderr, "E(x[%d]) = %g\n",
					i, objective_function(d, n, a, a[i]));

	float avg[d];
	linear_average(avg, d, n, a);
	E = objective_function(d, n, a, avg);
	fprintf(stderr, "avg = %lf %lf (%lf)\n", avg[0], avg[1], E);

	for (int j = 0; j < d; j++)
		x[j] = avg[j];
	fprintf(stderr, "starting energy = %g\n", E);

	int numit = NUMIT();
	for (int i = 0; i < numit; i++)
	{
		float g[d];
		gradient(g, d, n, a, x);
		if (d == 2)
			fprintf(stderr, "g[%d] = %g %g\t", i, g[0], g[1]);
		float λ = 0.5; // use weiszfeld weight here
		for (int j = 0; j < d; j++)
			x[j] -= λ * g[j];

		E = objective_function(d, n, a, x);
		fprintf(stderr, "f(%g %g) = %g\n", x[0], x[1], E);
	}

	return 0;
}

// fancier weiszfeld variands based on gradient descent
// (by default, plain weiszfeld)
//
// Variants:
// 	stepsize-constant: use a constant stepsize (relative to the gradient)
// 	stepsize-absolute: use a constant-length stepsize (non-convergent!)
// 	stepsize-factor: multiply the Weiszfeld lambda by this factor
// 	search-armijo: linear armijo search with parameters = 1/2
// 	stochastic-computations: 
//
int main_hess(int c, char *v[])
{
	char *out_sampling = pick_option(&c, &v, "o", "");

	if (c != 1)
		return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);

	int n, d;
	void *aa = iio_read_image_float("-", &d, &n);
	float (*a)[d] = aa;

	if (d == 2 && *out_sampling)
	{
		int w = 1000;
		int h = 1000;
		float *o = malloc(w*h*sizeof*o);
		for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			float x[2] = {
				XMIN() + i*(XMAX()-XMIN())/w,
				YMIN() + j*(YMAX()-YMIN())/h,
			};
			o[j*w+i] = objective_function(d, n, a, x);
		}
		iio_write_image_float(out_sampling, o, w, h);
		free(o);
	}

	float x[d]; // initialization
	for (int i = 0; i < d; i++)
		x[i] = 0;

	float E = objective_function(d, n, a, x);

	fprintf(stderr, "got %d points in dimension %d\n", n, d);
	fprintf(stderr, "energy at zero = %g\n", E);


	if (n < 10)
		for (int i = 0; i < n; i++)
			fprintf(stderr, "E(x[%d]) = %g\n",
					i, objective_function(d, n, a, a[i]));

	float avg[d];
	linear_average(avg, d, n, a);
	if (d == 2) {
		avg[0] = 9.99;
		avg[1] = 9.99;
	}
	E = objective_function(d, n, a, avg);
	fprintf(stderr, "avg = %lf %lf (%lf)\n", avg[0], avg[1], E);

	for (int j = 0; j < d; j++)
		x[j] = avg[j];
	fprintf(stderr, "starting energy = %g\n", E);

	int numit = NUMIT();
	for (int i = 0; i < numit; i++)
	{
		float g[d];
		gradient(g, d, n, a, x);
		if (d == 2)
			fprintf(stderr, "g[%d] = %g %g\n", i, g[0], g[1]);

		float H[d*d];
		hessian(H, d, n, a, x);
		if (d == 2) {
		fprintf(stderr, "\th=%g %g %g\n", H[0], H[1], H[3]);
		fprintf(stderr, "\tT=%g D=%g\n",H[0]+H[2], H[0]*H[3]-H[1]*H[2]);
		}

		float p[d];
		solve_spdf(p, H, g, d);
		if (d == 2) {
			fprintf(stderr, "\tp = %g %g\n", p[0], p[1]);
			float y[2] = {
				g[0] - H[0]*p[0] - H[1]*p[1],
				g[1] - H[2]*p[0] - H[3]*p[1],
			};
			fprintf(stderr, "\ty = %g %g\n", y[0], y[1]);
		}

		float W = weiszfeld_weight(d, n, a, x);
		fprintf(stderr, "\tW=%g, 1/W=%g\n", W, 1/W);


		float λ = 1; // use weiszfeld weight here
		for (int j = 0; j < d; j++)
			//x[j] -= λ * g[j];
			x[j] -= p[j];

		E = objective_function(d, n, a, x);
		fprintf(stderr, "\tf(%lf %lf) = %lf\n", x[0], x[1], E);

	}

	return 0;
}

#ifndef HIDE_ALL_MAINS
//int main(int c, char **v) { return main_clmps(c, v); }
//int main(int c, char **v) { return main_weisz(c, v); }
int main(int c, char **v) { return main_hess(c, v); }
#endif
