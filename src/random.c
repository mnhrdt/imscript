#ifndef _RANDOM_C
#define _RANDOM_C




#ifdef RNG_SPECIFIC_WELL1024
#include "well1024.c"
#else
#include <stdlib.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


static void xsrand(unsigned int seed)
{
#ifdef RNG_SPECIFIC_WELL1024
	well1024_seed(seed);
#else
#  ifndef __OpenBSD__
	srand(seed);
#  endif
#endif
}

static int xrand(void)
{
#ifdef RNG_SPECIFIC_WELL1024
	int r;
# ifdef _OPENMP
# pragma omp critical
# endif
	r = RAND_MAX * well1024();
	return r;
#else
#ifdef __OpenBSD__
	return arc4random();
#else
	return rand();
#endif
#endif
}

static double random_raw(void)
{
	return xrand();
}

static double random_uniform(void)
{
#ifdef RNG_SPECIFIC_WELL1024
	double r;
# ifdef _OPENMP
# pragma omp critical
# endif
	r = well1024();
	return r;
#else
#ifdef __OpenBSD__
	return arc4random_uniform(RAND_MAX)/(0.0+RAND_MAX);
#else
	return rand()/(1.0+RAND_MAX);
#endif
#endif
}

static double random_ramp(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	//double y = sqrt(random_uniform());
	double y = fmax(x1, x2);
	return y;
}

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y;
}

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
#ifdef __OpenBSD__
	return a + arc4random_uniform(b - a + 1);
#else
	return a + rand()%(b - a + 1);
#endif
}

static double random_laplace(void)
{
	double x = random_uniform();
	double y = random_uniform();
	double r = log(x/y);
	return isfinite(r)?r:0;
}

static double random_cauchy(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y_1 = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	double y_2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	double r = y_1/y_2;
	return isfinite(r)?r:0;
}

static double random_exponential(void)
{
	//double u = random_uniform();
	//double r = -log(1-u);
	//return r;
	return fabs(random_laplace());
}

static double random_pareto(void)
{
	return exp(random_exponential());
}

// This function samples a stable variable of parameters (alpha,beta)
//
// The algorithm is copied verbatim from formulas (2.3) and (2.4) in
//     Chambers, J. M.; Mallows, C. L.; Stuck, B. W. (1976).
//     "A Method for Simulating Stable Random Variables".
//     Journal of the American Statistical Association 71 (354): 340–344.
//     doi:10.1080/01621459.1976.10480344.
//
// Observation: the algorithm is numerically imprecise when alpha approaches 1.
// TODO: implement appropriate rearrangements as suggested in the article.
static double random_stable(double alpha, double beta)
{
	double U = (random_uniform() - 0.5) * M_PI;
	double W = random_exponential();
	double z = -beta * tan(M_PI * alpha / 2);
	double x = alpha == 1 ? M_PI/2 : atan(-z) / alpha;
	if (alpha == 0) {
		double a = (M_PI/2 + beta * U) * tan(U);
		double b = log(((M_PI/2) * W * cos(U)) / ((M_PI/2) + beta * U));
		return (a - beta * b) / x;
	} else {
		double a = pow(1 + z * z, 1 / (2*alpha));
		double b = sin(alpha * (U + x)) / pow(cos(U), 1/alpha);
		double c = pow(cos(U - alpha*(U + x)) / W, (1 - alpha) / alpha);
		return a * b * c;
	}
}

#endif//_RANDOM_C
