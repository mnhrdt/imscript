// implementation of Cohen-Lee-Miller-Pachocki-Sidford geometric medians
// ref. "Geometric Median in Nearly Linear Time", arXiv:160605225v

#include <assert.h>
#include <math.h>

#include "linalg.c" // cholesky, solve_spd


static double *global_log; // nine columns
static int global_i; // iteration counter for the log
static int global_maxit;
static double *global_solution;
static double global_optimal_energy;
#define LOG_i 0
#define LOG_xi 1
#define LOG_yi 2
#define LOG_fi 3
#define LOG_gxi 4
#define LOG_gyi 5
#define LOG_gni 6
#define LOG_di 7
#define LOG_ei 8

static void linear_average(double *o, int d, int n, double a[n][d])
{
	for (int k = 0; k < d; k++)
	{
		o[k] = 0;
		for (int i = 0; i < n; i++)
			o[k] += a[i][k] / n;
	}
}

// y[k] = (1/n) * sum_i x[i][k]
static void double_avg(double *y, double *xx, int d, int n)
{
	double (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

static double fnorm(double *x, int n)
{
	switch(n) {
	case 1: return fabs(x[0]);
	case 2: return hypot(x[0], x[1]);
	default: return hypot(x[0], fnorm(x+1, n-1));
	}
}

static double scalar_product(int n, double x[n], double y[n])
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

static double medscore(double *xx, int idx, int d, int n)
{
	double (*x)[d] = (void*)xx;
	double r = 0;
	for (int i = 0; i < n; i++)
		if (i != idx)
		{
			double v[d];
			for (int j = 0; j < d; j++)
				v[j] = x[idx][j] - x[i][j];
			r += fnorm(v, d);
		}
	return r;
}

// y[] = x[i][] which is closest to the euclidean median
// runs in quadratic time
static void double_med(double *y, double *xx, int d, int n)
{
	double (*x)[d] = (void*)xx;
	int midx = 0;
	double misc = medscore(xx, 0, d, n);
	for (int i = 1; i < n; i++)
	{
		double si = medscore(xx, i, d, n);
		if (si < misc) {
			midx = i;
			misc = si;
		}
	}
	for (int i = 0; i < d; i++)
		y[i] = x[midx][i];

}

static double euclidean_distance(double *p, double *q, int d)
{
	double r = 0;
	for (int k = 0; k < d; k++)
		r = hypot(r, p[k] - q[k]);
	return r;
}

static double objective_function(int d, int n, double a[n][d], double *x)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += euclidean_distance(a[i], x, d);
	return r;
}

static double penalized_objective(double *a[], int d, int n, double *x,
		double e, double t)
{
	double r = 0;
	for (int i = 0; i < n; i++)
	{
		double e = euclidean_distance(x, a[i], d);
		double s = sqrt(1 + t*t * e*e);
		r += s - log(1 + s);
	}
	return r;
}

void clmps_median(
		double *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		double a[n][d],    // input list of points
		double e           // desired accuracy 0<e<1
		)
{
	assert(0 < e && e < 1);

	// compute a 2-approximate geometric median and use it to center
	double x[d];
	linear_average(x, d, n, a);
	double f_star = objective_function(d, n, a, x);
	double t[n];
	for (int i = 0; i < n; i++)
		t[i] = pow(601.0/600, i)/(400*f_star);
	double e_star = e/3;
	double t_star = 2 * n / (e_star * f_star);
	double e_v = pow(e_star/(7*n), 2)/8;
	double e_c = pow(e_v/36, 1.5);
	double y[d];
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

// XXX warning! non-reentrant
static void reservoir_sampling_non_reentrant(int *s, int n, int k)
{
	// "Algorithm L", Kim-Hung Li, 1994
	// O(k + k*log(n/k))
	static int *r = NULL;
	static int rn = 0;
	if (n != rn) {
		if (r) free(r);
		if (n) r = malloc(n * sizeof*r);
		rn = n;
	}

	for (int i = 0; i < k; i++)
		r[i] = s[i];

	double W = exp(log(random_uniform())/k);

	int i = k;
	while (i < n)
	{
		i = i + floor(log(random_uniform())/log(1-W));
		if (i < n) {
			r[randombounds(0,k-1)] = s[i];
			W = W*exp(log(random_uniform())/k);
		}
	}
}

static int compare_doubles(const void *aa, const void *bb)
{
	const double *a = (const double *)aa;
	const double *b = (const double *)bb;
	return (*a > *b) - (*a < *b);
}

double crude_approximate(
		double *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		double a[n][d],    // input list of points
		int K             // subset size
		)
{
	assert(K > 0 && K < n);

	int S1[K], S2[K];
	get_random_subset(S1, n, K);
	get_random_subset(S2, n, K);

	double α[K][K];
	for (int i = 0; i < K; i++)
	for (int j = 0; j < K; j++)
		α[i][j] = euclidean_distance(a[S2[i]], a[S1[j]], d);

	for (int i = 0; i < K; i++)
		qsort(α[i], K, sizeof*α[i], compare_doubles);

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
		double *out,       // output median (vector of length d)
		int d,            // dimension of each point
		int n,            // number of input points
		double a[n][d],    // input list of points
		double ε           // desired precision
		)
{
	double T = pow(60/ε, 2);
	double x[d];
	double λ = crude_approximate(x, d, n, a, lrint(sqrt(T)));
	double η = 6*λ*sqrt(2/T)/n;
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
		double ρ = euclidean_distance(x, a[i], d);
		double g[d];
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
static double fdiste(double *x, double *y, int n, double e)
{
	return n ? hypot(*x - *y, fdiste(x + 1, y + 1, n - 1, e)) : e;
}

//#include "smapa.h"
//SMART_PARAMETER_SILENT(WEISZ_NITER,10)

// y[k] = euclidean median of the vectors x[i][k]
static void double_weisz(double *y, double *x, int d, int n, int N, double ε)
{
	double_avg(y, x, d, n);
	for (int k = 0; k < N; k++) {
		double a[d], b = 0;
		for (int l = 0; l < d; l++)
			a[l] = 0;
		for (int i = 0; i < n; i++) {
			double dxy = fdiste(x + i*d, y, d, ε);
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

static void fill_points(int d, int n, double x[n][d])
{
	assert(d == 2);

	int m = 3;         // number of centers
	double c[][2] = {  // centers
		{100, 240},
		{400, 250},
		{100, 260},
	};
	double σ = 4;

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

static void fill_points2(int d, int n, double x[n][d])
{
	assert(d == 2);

	int m = 2;         // number of centers
	double c[][2] = {  // centers
		{100, 100},
		{400, 400}
	};
	double σ[] = {15, 3};

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
	double weisz_n = atoi(pick_option(&c, &v, "wn", "10"));
	double weisz_e = atof(pick_option(&c, &v, "we", "0.01"));
	double crude_K = atoi(pick_option(&c, &v, "ck", "10"));
	double linap_e = atof(pick_option(&c, &v, "le", "0.5"));

	int d = 2;     // dimension
	int n = 3*10000; // total number of points
	double x[n][d];
	fill_points2(d, n, x);

	double E;

	double avg[d];
	linear_average(avg, d, n, x);
	E = objective_function(d, n, x, avg)/n;
	fprintf(stderr, "avg = %lf %lf (%lf)\n", avg[0], avg[1], E);

	double wei[d];
	double_weisz(wei, *x, d, n, weisz_n, weisz_e);
	E = objective_function(d, n, x, wei)/n;
	fprintf(stderr, "wei = %lf %lf (%lf)\n", wei[0], wei[1], E);

	double cru[d];
	double r = crude_approximate(cru, d, n, x, crude_K);
	E = objective_function(d, n, x, cru)/n;
	fprintf(stderr, "cru = %lf %lf {%lf} (%lf)\n", cru[0], cru[1], r, E);

	double lin[d];
	linear_approximate_median(lin, d, n, x, linap_e);
	E = objective_function(d, n, x, lin)/n;
	fprintf(stderr, "lin = %lf %lf (%lf)\n", lin[0], lin[1], E);

	double med[d];
	double_med(med, *x, d, n);
	E = objective_function(d, n, x, med)/n;
	fprintf(stderr, "med = %lf %lf (%lf)\n", med[0], med[1], E);


	return 0;
}

//static double objective_function(int d, int n, double a[n][d], double *x)
//{
//	double r = 0;
//	for (int i = 0; i < n; i++)
//		r += euclidean_distance(a[i], x, d);
//	return r;
//}

static void gradient(double *g, int d, int n, double a[n][d], double *x)
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

static double weiszfeld_weight(int d, int n, double a[n][d], double *x)
{
	double R = 0; // sum of all inverse norms x-ai
	for (int i = 0; i < n; i++)
	{
		double r = 0; // norm of x-ai
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);
		R += 1/r;
	}
	return 1/R;
}

static void hessian(double *h, int d, int n, double a[n][d], double *x)
{
	double (*H)[d] = (void*)h;

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

	//if (d == 2)
	//	fprintf(stderr, "\t\tH=%g %g %g %g\n", h[0],h[1],h[2],h[3]);

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
}

// compute only the diagonal of the hessian
static void hessian_diag(double *h, int d, int n, double a[n][d], double *x)
{
	// initialize to the identity
	for (int k = 0; k < d; k++)
		h[k] = k;

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
		h[k] *= R;

	// second-term: accumulate singular perturbations
	for (int i = 0; i < n; i++)
	{
		double r = 0; // norm of x-ai (TODO: do not recompute)
		for (int j = 0; j < d; j++)
			r = hypot(r, x[j] - a[i][j]);

		for (int k = 0; k < d; k++)
			h[k] -= (x[k]-a[i][k])*(x[k]-a[i][k]) / (r*r*r);
	}
}

#include "iio.h"

//#include "smapa.h"
//SMART_PARAMETER_SILENT(XMIN,-2)
//SMART_PARAMETER_SILENT(XMAX,12)
//SMART_PARAMETER_SILENT(YMIN,-2)
//SMART_PARAMETER_SILENT(YMAX,12)
//SMART_PARAMETER_SILENT(NUMIT,20)

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
//int main_weisz(int c, char *v[])
//{
//	char *out_sampling = pick_option(&c, &v, "o", "");
//
//	if (c != 1)
//		return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);
//
//	int n, d;
//	void *aa = iio_read_image_double("-", &d, &n);
//	double (*a)[d] = aa;
//
//	if (d == 2 && *out_sampling)
//	{
//		int w = 800;
//		int h = 600;
//		double *o = malloc(w*h*sizeof*o);
//		for (int i = 0; i < h; i++)
//		for (int j = 0; j < w; j++)
//		{
//			double x[2] = {
//				XMIN() + i*(XMAX()-XMIN())/w,
//				YMIN() + j*(YMAX()-YMIN())/h,
//			};
//			o[j*w+i] = objective_function(d, n, a, x);
//		}
//		iio_write_image_double(out_sampling, o, w, h);
//		free(o);
//	}
//
//	double x[d]; // initialization
//	for (int i = 0; i < d; i++)
//		x[i] = 0;
//
//	double E = objective_function(d, n, a, x);
//
//	fprintf(stderr, "got %d points in dimension %d\n", n, d);
//	fprintf(stderr, "energy at zero = %g\n", E);
//
//
//	if (n < 10)
//		for (int i = 0; i < n; i++)
//			fprintf(stderr, "E(x[%d]) = %g\n",
//					i, objective_function(d, n, a, a[i]));
//
//	double avg[d];
//	linear_average(avg, d, n, a);
//	E = objective_function(d, n, a, avg);
//	fprintf(stderr, "avg = %lf %lf (%lf)\n", avg[0], avg[1], E);
//
//	for (int j = 0; j < d; j++)
//		x[j] = avg[j];
//	fprintf(stderr, "starting energy = %g\n", E);
//
//	int numit = NUMIT();
//	for (int i = 0; i < numit; i++)
//	{
//		double g[d];
//		gradient(g, d, n, a, x);
//		if (d == 2)
//			fprintf(stderr, "g[%d] = %g %g\t", i, g[0], g[1]);
//		double λ = 0.5; // use weiszfeld weight here
//		for (int j = 0; j < d; j++)
//			x[j] -= λ * g[j];
//
//		E = objective_function(d, n, a, x);
//		fprintf(stderr, "f(%g %g) = %g\n", x[0], x[1], E);
//	}
//
//	return 0;
//}

// halve the step until we first observe an improvement
static void line_search_until_first_decrease(
		double *y,        // output position
		int d,           // dimension
		int n,           // number of data points
		double a[n][d],   // input data points
		double *x,        // base point position
		double *p,        // minus descent direction
		double t          // initial overshot
		)
{
	assert(t > 0);
	double f0 = objective_function(d, n, a, x);
	while (t > 1e-6)
	{
		for (int k = 0; k < d; k++)
			y[k] = x[k] - t * p[k];
		double ft = objective_function(d, n, a, y);
		if (ft < f0)
			break;
		t *= 0.5; // armijo constant
	}
	fprintf(stderr, "\tlsufd t = %g\n", t);
}

// backtrack until first decrease
// TODO: make it fill-in the new point instead of returning the step
static double backtrack_first(
		int d,           // dimension
		int n,           // number of data points
		double a[n][d],  // input data points
		double *x,       // base point position
		double *p,       // minus descent direction
		double t         // initial overshot
		)
{
	assert(t > 0);
	double f0 = objective_function(d, n, a, x);
	while (t > 1e-7)
	{
		double y[d];
		for (int k = 0; k < d; k++)
			y[k] = x[k] - t * p[k];
		double ft = objective_function(d, n, a, y);
		if (ft < f0)
			break;
		t *= 0.5; // armijo constant
	}
	return t;
}

// halve the step until armijo condition holds
static double backtrack_armijo(
		int d,           // dimension
		int n,           // number of data points
		double a[n][d],   // input data points
		double *x,        // base point position
		double *p,        // minus descent direction
		double t          // initial overshot
		)
{
	assert(t > 0);
	double f0 = objective_function(d, n, a, x);
	double g[d];
	gradient(g, d, n, a, x);
	double m = -scalar_product(d, g, p);
	assert(m < 0);
	double σ = 0.5; // second Armijo constant

	while (t > 1e-6)
	{
		double y[d];
		for (int k = 0; k < d; k++)
			y[k] = x[k] - t * p[k];
		double ft = objective_function(d, n, a, y);
		if (ft < f0 + σ * t * m)
			break;
		t *= 0.5; // first Armijo constant
	}
	return t;
}

// forward search until last decrease
static double forward_last(
		int d,           // dimension
		int n,           // number of data points
		double a[n][d],  // input data points
		double *x,       // base point position
		double *p,       // minus descent direction
		double t         // initial undershot
		)
{
	assert(t > 0);
	double inverse_armijo_constant = 2;
	double f_old = objective_function(d, n, a, x);
	double t_old = t;
	while (1)
	{
		double y[d];
		for (int k = 0; k < d; k++)
			y[k] = x[k] - t * p[k];
		double f = objective_function(d, n, a, y);
		if (f >= f_old)
			return t_old;
		else {
			t_old = t;
			f_old = f;
			t *= inverse_armijo_constant;
		}
	}
	return t;
}

// halve the step until armijo condition holds
static void line_search_armijo(
		double *y,        // output position
		int d,           // dimension
		int n,           // number of data points
		double a[n][d],   // input data points
		double *x,        // base point position
		double *p,        // minus descent direction
		double t          // initial overshot
		)
{
	assert(t > 0);
	double f0 = objective_function(d, n, a, x);
	double g[d];
	gradient(g, d, n, a, x);
	double m = -scalar_product(d, g, p);
	assert(m < 0);
	double σ = 0.5; // second Armijo constant

	while (t > 1e-6)
	{
		for (int k = 0; k < d; k++)
			y[k] = x[k] - t * p[k];
		double ft = objective_function(d, n, a, y);
		if (ft < f0 + σ * t * m)
			break;
		t *= 0.5; // first A/rmijo constant
	}
	fprintf(stderr, "\tlsufd t = %g\n", t);
}

struct descent_options {

	// -i
	enum {
		ZERO,             // initialize at x_0=0
		GIVEN,            // user-provided x_0
		AVERAGE,          // x_0 = avg(a_1, ... a_n)
		GORNER_KANZOW     // argmin(f(a_i)) - λ_k grad(f)  (very slow!)
	} initialization;

	// -d
	enum {
		GRADIENT,         // gradient
		NEWTON,           // hessian \ gradient
		DIAGONAL_HESSIAN  // diagonal of the hessian \ gradient
	} descent_direction;

	// -s
	enum {
		RELATIVE,         // raw descent direction
		ABSOLUTE,         // normalize descent direction
		WEISZFELD,        // weiszfeld scaling factor
		BACK_FIRST,       // backtrack until first objective lowering
		BACK_ARMIJO,      // backtrack until armijo condition
		FWD_LAST,         // advance while f decreases
		FWD_ARMIJO,       // advance until armijo
		BFGS              // BFGS
	} descent_step;

	// (defined by -p)
	enum {
		FULL,              // use whole dataset
		PROPORTION,        // use a fraction of the dataset
		AMOUNT             // use a fixed amount of data points
	} stochasticity; // note: points chosen using reservoir sampling

	// -l
	double step_size;

	// -p
	double stochasticity_parameter;

	//double step_factor;            // 0.5
	//double second_armijo_constant; // 0.5
};

static void print_descent_options(FILE *f, struct descent_options *o)
{
	char const*const s_i[] = {
		[ZERO] = "ZERO",
		[GIVEN] = "GIVEN",
		[AVERAGE] = "AVERAGE",
		[GORNER_KANZOW] = "GORNER_KANZOW",
	};
	char const*const s_d[] = {
		[GRADIENT] = "GRADIENT",
		[NEWTON] = "NEWTON",
		[DIAGONAL_HESSIAN] = "DIAGONAL_HESSIAN",
	};
	char const*const s_s[] = {
		[RELATIVE] = "RELATIVE",
		[ABSOLUTE] = "ABSOLUTE",
		[WEISZFELD] = "WEISZFELD",
		[BACK_FIRST] = "BACK_FIRST",
		[BACK_ARMIJO] = "BACK_ARMIJO",
		[FWD_LAST] = "FWD_LAST",
		[FWD_ARMIJO] = "FWD_ARMIJO",
		[BFGS] = "BFGS",
	};
	char const*const s_p[] = {
		[FULL] = "FULL",
		[PROPORTION] = "PROPORTION",
		[AMOUNT] = "AMOUNT",
	};
	fprintf(f, "O initialization = %s\n", s_i[o->initialization]);
	fprintf(f, "O direction = %s\n", s_d[o->descent_direction]);
	fprintf(f, "O step = %s\n", s_s[o->descent_step]);
	fprintf(f, "O stochasticity = %s\n", s_p[o->stochasticity]);
	fprintf(f, "O step_size = %g\n", o->step_size);
	fprintf(f, "O stoch_parameter = %g\n", o->stochasticity_parameter);
}

static void find_initial_point(
		int d,            // dimension of the space
		int n,            // number of input points
		double a[n][d],   // input array of points
		double x[d],      // output initial point
		struct descent_options *o
		)
{
	switch(o->initialization) {
	case ZERO:
		for (int i = 0; i < d; i++)
			x[i] = 0;
		break;
	case AVERAGE:
		linear_average(x, d, n, a);
		break;
	case GORNER_KANZOW:
		exit(fprintf(stderr, "Gorner Kanzow init not implemented\n"));
		break;
	default:
		exit(fprintf(stderr, "bad init (%d)\n", o->initialization));
		break;
	}
}

static void find_descent_direction(
		int d,            // dimension of the space
		int n,            // number of input points
		double p[d],      // output descent direction
		double a[n][d],   // input array of points
		double x[d],      // input base point
		struct descent_options *o
		)
{
	switch(o->descent_direction) {
	case GRADIENT:
		gradient(p, d, n, a, x);
		break;
	case NEWTON:
		{
			double g[d], H[d*d];
			gradient(g, d, n, a, x);
			hessian(H, d, n, a, x);
			solve_spd(p, H, g, d);
			break;
		}
	case DIAGONAL_HESSIAN:
		{
			double g[d], h[d];
			gradient(g, d, n, a, x);
			hessian_diag(h, d, n, a, x);
			for (int i = 0; i < d; i++)
				p[i] = g[i] / h[i];
			break;
		}
	}
}

// find the factor by which to scale p to build the next x
static double find_descent_step(
		int d,            // dimension of the space
		int n,            // number of input points
		double p[d],      // input descent direction
		double a[n][d],   // input array of points
		double x[d],      // input base point
		struct descent_options *o
		)
{
	// TODO: add BFGS criterion
	// see e.g. stoer-bulirsch chapter on unconstrained optimization
	switch(o->descent_step) {
	case RELATIVE:    return o->step_size;
	case ABSOLUTE:    return o->step_size / fnorm(p, d);
	case WEISZFELD:   return o->step_size * weiszfeld_weight(d, n, a, x);
	case BACK_FIRST:  return backtrack_first (d, n, a, x, p, o->step_size);
	case BACK_ARMIJO: return backtrack_armijo(d, n, a, x, p, o->step_size);
	case FWD_LAST:    return forward_last    (d, n, a, x, p, o->step_size);
	//case FWD_ARMIJO:  return forward_armijo  (d, n, a, x, p, o->step_size);
	default:          return 1;
	}
}

static void grab_options(struct descent_options *o, int *c, char ***v)
{
	char *opt_i = pick_option(c, v, "i", "average");
	char *opt_d = pick_option(c, v, "d", "gradient");
	char *opt_s = pick_option(c, v, "s", "weiszfeld");
	char *opt_l = pick_option(c, v, "l", "1.0");
	char *opt_p = pick_option(c, v, "p", "0");

	if (!strcmp(opt_i, "zero"))    o->initialization = ZERO;
	if (!strcmp(opt_i, "average")) o->initialization = AVERAGE;
	if (!strcmp(opt_i, "gk"))      o->initialization = GORNER_KANZOW;
	//if (!strcmp(opt_i, "given"))   o->initialization = GIVEN;
	if (!strcmp(opt_d, "gradient")) o->descent_direction = GRADIENT;
	if (!strcmp(opt_d, "newton"))   o->descent_direction = NEWTON;
	if (!strcmp(opt_d, "dhessian")) o->descent_direction = DIAGONAL_HESSIAN;
	if (!strcmp(opt_s, "relative"))    o->descent_step = RELATIVE;
	if (!strcmp(opt_s, "absolute"))    o->descent_step = ABSOLUTE;
	if (!strcmp(opt_s, "weiszfeld"))   o->descent_step = WEISZFELD;
	if (!strcmp(opt_s, "back-first"))  o->descent_step = BACK_FIRST;
	if (!strcmp(opt_s, "back-armijo")) o->descent_step = BACK_ARMIJO;
	if (!strcmp(opt_s, "fwd-last"))    o->descent_step = FWD_LAST;
	if (!strcmp(opt_s, "fwd-armijo"))  o->descent_step = FWD_ARMIJO;

	o->step_size = atof(opt_l);

	// parse the stochasticity parameters
	int n = strlen(opt_p);
	double p = atof(opt_p);
	o->stochasticity = FULL;
	if (p > 0 && p < 1) {
		o->stochasticity = PROPORTION;
		o->stochasticity_parameter = p;
	}
	if (p > 0 && p <= 100 && opt_p[n-1] == '%') {
		o->stochasticity = PROPORTION;
		o->stochasticity_parameter = p/100;
	}
	if (p > 0 && opt_p[n-1] != '%') {
		o->stochasticity = AMOUNT;
		o->stochasticity_parameter = p;
	}
	//exit(fprintf(stderr, "p=%g percentage=%d\n", p, opt_p[n-1]=='%'));
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
//int main_hess(int c, char *v[])
//{
//	struct descent_options o[1];
//	grab_options(o, &c, &v);
//	char *out_sampling = pick_option(&c, &v, "o", "");
//	char *out_log = pick_option(&c, &v, "v", "");
//	//char *known_solution_s = pick_option(&v, &v, "s", "");
//
//	if (c != 1)
//		return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);
//
//	int n, d;
//	void *aa = iio_read_image_double("-", &d, &n);
//	double (*a)[d] = aa;
//
//	if (d == 2 && *out_sampling)
//	{
//		int w = 800;
//		int h = 600;
//		double *o = malloc(w*h*sizeof*o);
//		for (int j = 0; j < h; j++)
//		for (int i = 0; i < w; i++)
//		{
//			double x[2] = {
//				XMIN() + i*(XMAX()-XMIN())/w,
//				YMIN() + j*(YMAX()-YMIN())/h,
//			};
//			o[j*w+i] = objective_function(d, n, a, x);
//		}
//		iio_write_image_double(out_sampling, o, w, h);
//		free(o);
//	}
//
//	if (*out_log)
//	{
//		global_i = 0;
//		global_maxit = NUMIT();
//		global_log = malloc(9 * global_maxit * sizeof*global_log);
//		global_solution = malloc(d * sizeof*global_solution);
//		for (int i = 0; i < d; i++)
//			global_solution[i] = NAN;
//	}
//
//	double x[d]; // initialization
//	for (int i = 0; i < d; i++)
//		x[i] = 0;
//
//	double E = objective_function(d, n, a, x);
//
//	fprintf(stderr, "got %d points in dimension %d\n", n, d);
//	fprintf(stderr, "energy at zero = %g\n", E);
//
//
//	if (n < 10)
//		for (int i = 0; i < n; i++)
//			fprintf(stderr, "E(x[%d]) = %g\n",
//					i, objective_function(d, n, a, a[i]));
//
//	double avg[d];
//	linear_average(avg, d, n, a);
//	if (d == 2) {
//		//avg[0] = 9.99;
//		//avg[1] = 9.99;
//		avg[0] = 8;
//		avg[1] = 0;
//	}
//	E = objective_function(d, n, a, avg);
//	fprintf(stderr, "start = %lf %lf (%lf)\n", avg[0], avg[1], E);
//
//	for (int j = 0; j < d; j++)
//		x[j] = avg[j];
//	fprintf(stderr, "starting energy = %g\n", E);
//
//	int numit = NUMIT();
//	for (int i = 0; i < numit; i++)
//	{
//		double g[d];
//		gradient(g, d, n, a, x);
//		if (d == 2)
//			fprintf(stderr, "g[%d] = %g %g\n", i, g[0], g[1]);
//
//		double H[d*d];
//		hessian(H, d, n, a, x);
//		if (d == 2) {
//		fprintf(stderr, "\th=%g %g %g\n", H[0], H[1], H[3]);
//		fprintf(stderr, "\tT=%g D=%g\n",H[0]+H[2], H[0]*H[3]-H[1]*H[2]);
//		}
//
//		double p[d];
//		solve_spd(p, H, g, d);
//		if (d == 2) {
//			fprintf(stderr, "\tp = %g %g\n", p[0], p[1]);
//			double y[2] = {
//				g[0] - H[0]*p[0] - H[1]*p[1],
//				g[1] - H[2]*p[0] - H[3]*p[1],
//			};
//			fprintf(stderr, "\ty = %g %g\n", y[0], y[1]);
//		}
//
//		double W = weiszfeld_weight(d, n, a, x);
//		fprintf(stderr, "\tW=%g, 1/W=%g\n", W, 1/W);
//
//		double y[d];
//		line_search_until_first_decrease(y, d, n, a, x, p, 1.0);
//		//line_search_armijo(y, d, n, a, x, p, 1.0);
//		for (int k = 0; k < d; k++) x[k] = y[k];
//
//		//double λ = W; // use weiszfeld weight here
//		//for (int j = 0; j < d; j++)
//		//	x[j] -= λ * g[j];
//		//	x[j] -= p[j];
//
//
//		E = objective_function(d, n, a, x);
//		fprintf(stderr, "\tf(%lf %lf) = %lf\n", x[0], x[1], E);
//
//		if (*out_log && d > 1)
//		{
//			global_i = i;
//			global_log[9*global_i + LOG_i  ] = global_i;
//			global_log[9*global_i + LOG_xi ] = x[0];
//			global_log[9*global_i + LOG_yi ] = x[1];
//			global_log[9*global_i + LOG_fi ] = E;
//			global_log[9*global_i + LOG_gxi] = g[0];
//			global_log[9*global_i + LOG_gxi] = g[1];
//			global_log[9*global_i + LOG_gni] = hypot(g[0],g[1]);
//			global_log[9*global_i + LOG_ei ] = NAN;;
//			global_log[9*global_i + LOG_di ] = NAN;;
//		}
//
//	}
//
//	if (*out_log)
//	{
//		FILE *f = fopen(out_log, "w");
//		for (int i = 0; i < global_i; i++)
//		for (int j = 0; j < 9; j++)
//			fprintf(f, "%lf%c", global_log[9*i+j], j==8?'\n':' ');
//		fclose(f);
//	}
//
//	return 0;
//}


static char *help_string_name     = "geomedian";
static char *help_string_version  = "geomedian 1.0\n\nWritten by eml";
static char *help_string_oneliner = "algorithms for the geometric median";
static char *help_string_usage    = "usage:\n\t"
"geomedian [options] < points.txt > median.txt";
static char *help_string_long     =
"Geomedian approximates the geometric median of a point cloud.\n"
"\n"
"Various algorithms can be used depending on the options.\n"
"\n"
"Usage: geomedian [options] < points.txt > median.txt\n"
"\n"
"Options:\n"
" -n n\t\ttotal number of descent iterations\n"
" -i {zero,average,gk}\t\tselect initial point\n"
" -d {gradient,newton,dhessian}\t\tdescent direction\n"
"\n"
"Report bugs to <enric.meinhardt@ens-paris-saclay.fr>."
;
#include "help_stuff.c" // functions that print the strings named above

// general descent method, with weiszfeld and newton-armijo as particular cases
int main_descent(int c, char *v[])
{
	// process "help" arguments
	if (c == 2)
		if_help_is_requested_print_it_and_exit_the_program(v[1]);

	// extract named arguments
	struct descent_options o[1];
	grab_options(o, &c, &v);
	int numit = atoi(pick_option(&c, &v, "n", "100"));
	//char *out_log = pick_option(&c, &v, "v", "");      // for debugging

	// no positional arguments
	if (c != 1)
		return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);

	// read input point cloud
	int d; // dimension of the space
	int n; // total number of input points
	void *aa = iio_read_image_double("-", &d, &n);
	double (*a)[d] = aa; // point cloud

	//// debugging stuff
	//if (*out_log)
	//{
	//	global_i = 0;
	//	global_maxit = NUMIT();
	//	global_log = malloc(9 * global_maxit * sizeof*global_log);
	//	global_solution = malloc(d * sizeof*global_solution);
	//	for (int i = 0; i < d; i++)
	//		global_solution[i] = NAN;
	//}

	double x[d]; // current position
	find_initial_point(d, n, a, x, o);

	double E = objective_function(d, n, a, x);

	print_descent_options(stderr, o);
	fprintf(stderr, "got %d points in dimension %d\n", n, d);
	fprintf(stderr, "energy at first point = %g\n", E);
	fprintf(stderr, "numit = %d\n", numit);

	if (d == 2)
		fprintf(stderr, "P\t%lf\t%lf\t%lf\n", x[0], x[1], E);

	//if (n < 10)
	//	for (int i = 0; i < n; i++)
	//		fprintf(stderr, "E(x[%d]) = %g\n",
	//				i, objective_function(d, n, a, a[i]));

	for (int i = 0; i < numit; i++)
	{
		double p[d];
		find_descent_direction(d, n, p, a, x, o);

		double λ = find_descent_step(d, n, p, a, x, o);

		for (int k = 0; k < d; k++)
			x[k] = x[k]  -  λ * p[k];

		E = objective_function(d, n, a, x);
		if (d == 2)
			fprintf(stderr, "P\t%lf\t%lf\t%lf\n", x[0], x[1], E);

		//if (*out_log && d > 1)
		//{
		//	global_i = i;
		//	global_log[9*global_i + LOG_i  ] = global_i;
		//	global_log[9*global_i + LOG_xi ] = x[0];
		//	global_log[9*global_i + LOG_yi ] = x[1];
		//	global_log[9*global_i + LOG_fi ] = E;
		//	global_log[9*global_i + LOG_gxi] = g[0];
		//	global_log[9*global_i + LOG_gxi] = g[1];
		//	global_log[9*global_i + LOG_gni] = hypot(g[0],g[1]);
		//	global_log[9*global_i + LOG_ei ] = NAN;;
		//	global_log[9*global_i + LOG_di ] = NAN;;
		//}

	}

	//if (*out_log)
	//{
	//	FILE *f = fopen(out_log, "w");
	//	for (int i = 0; i < global_i; i++)
	//	for (int j = 0; j < 9; j++)
	//		fprintf(f, "%lf%c", global_log[9*i+j], j==8?'\n':' ');
	//	fclose(f);
	//}

	return 0;
}

int main_descent_raw(int c, char *v[])
{
	struct descent_options o[1];
	grab_options(o, &c, &v);
	int numit = atoi(pick_option(&c, &v, "n", "100"));
	if (c != 1) return fprintf(stderr, "usage:\n\t%s [params] <in\n", *v);

	// read input point cloud
	int d; // dimension of the space
	int n; // total number of input points
	void *aa = iio_read_image_double("-", &d, &n);
	double (*a)[d] = aa; // point cloud

	double x[d]; // current position
	find_initial_point(d, n, a, x, o);

	double E = objective_function(d, n, a, x);

	for (int i = 0; i < numit; i++)
	{
		double p[d];
		find_descent_direction(d, n, p, a, x, o);

		double λ = find_descent_step(d, n, p, a, x, o);

		for (int k = 0; k < d; k++)
			x[k] = x[k]  -  λ * p[k];

	//	E = objective_function(d, n, a, x);
	}
	return 0;
}

#ifndef HIDE_ALL_MAINS
//int main(int c, char **v) { return main_clmps(c, v); }
//int main(int c, char **v) { return main_weisz(c, v); }
//int main(int c, char **v) { return main_hess(c, v); }
int main(int c, char **v) { return main_descent(c, v); }
#endif
