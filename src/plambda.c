// NAME
// 	plambda - a RPN calculator for image pixels
//
// SYNOPSIS
// 	plambda a.pnb b.png c.png ... "lambda expression" > output.png
//
// DESCRIPTION
// 	Plambda applies an expression to all the pixels of a collection of
// 	images, and produces a single output image.  Each input image
// 	corresponds to one of the variables of the expression (in alphabetical
// 	order).  There are modifiers to the variables that allow access to the
// 	values of neighboring pixels, or to particular components of a pixel.
//
// LANGUAGE
// 	A "plambda" program is a sequence of tokens.  Tokens may be constants,
// 	variables, or operators.  Constants and variables get their value
// 	computed and pushed to the stack.  Operators pop values from the stack,
// 	apply a function to them, and push back the results.
//
// 	CONSTANTS: numeric constants written in scientific notation, and "pi"
// 	OPERATORS: +, -, *, ^, /, and all the functions from math.h
// 	VARIABLES: anything not recognized as a constant or operator.  There
// 	must be as many variables as input images, and they are assigned to
// 	images in alphabetical order.
//
//	Some "sugar" is added to the language:
//
//	Predefined variables (always preceeded by a colon):
//
//		TOKEN	MEANING
//
//		:i	horizontal coordinate of the pixel
//		:j	vertical coordinate of the pixel
//		:w	width of the image
//		:h	heigth of the image
//		:n	number of pixels in the image
//		:x	relative horizontal coordinate of the pixel
//		:y	relative horizontal coordinate of the pixel
//		:r	relative distance to the center of the image
//		:t	relative angle from the center of the image
//		:I	horizontal coordinate of the pixel (centered)
//		:J	vertical coordinate of the pixel (centered)
//		:W	width of the image divided by 2*pi
//		:H	height of the image divided by 2*pi
//
//
//	Variable modifiers acting on regular variables:
//
//		TOKEN	MEANING
//
//		x	value of pixel (i,j)
//		x(0,0)	value of pixel (i,j)
//		x(1,0)	value of pixel (i+1,j)
//		x(0,-1)	value of pixel (i,j-1)
//		...
//
//		x	value of pixel (i,j)
//		x[0]	value of first component of pixel (i,j)
//		x[1]	value of second component of pixel (i,j)
//
//		x(1,-1)[2] value of third component of pixel (i+1,j-1)
//
//	Stack operators (allow direct manipulation of the stack):
//
//		TOKEN	MEANING
//
//		del	remove the value at the top of the stack (ATTOTS)
//		dup	duplicate the value ATTOTS
//		rot	swap the two values ATTOTS
//		split	split the vector ATTOTS into scalar components
//		join	join the components of two vectors ATTOTS
//		join3	join the components of three vectors ATTOTS
//
//	Magic modifiers (which can not be computed locally for each image):
//
//		x%i	value of the smallest sample of image "x"
//		x%a	value of the largest sample
//		x%v	average sample value
//		x%m	median sample value
//		x%I	value of the smallest pixel
//		x%A	value of the largest pixel
//		x%V	average pixel value
//		x%M	median pixel value (not implemented)
//		x%qn	nth sample percentile
//		x%Qn	nth pixel percentile (not implemented)
//		x%r	random sample of the image (not implemented)
//		x%R	random pixel of the image (not implemented)
//
//	Notice that the scalar "magic" modifiers may act upon individual
//	components, e.g. x[2]%i is the minimum value of the blue component.
//
//
//	Other operators:
//
//		TOKEN	MEANING
//		vmprod3x3 
//
// EXAMPLES
// 	Sum two images:
//
// 		plambda a.png b.png "a b +" > aplusb.png
//
//	Add a gaussian to half of lena:
//
//		plambda /tmp/lena.png "x 2 / :r :r * -1 * 40 * exp 200 * +"
//
//	Forward differences to compute the derivative in horizontal direction:
//
//		plambda lena.png "x(1,0) x -"
//
//	Sobel edge detector:
//		plambda lena.png "x(1,0) 2 * x(1,1) x(1,-1) + + x(-1,0) 2 * x(-1,1) x(-1,-1) + + - x(0,1) 2 * x(1,1) x(-1,1) + + x(0,-1) 2 * x(1,-1) x(-1,-1) + + - hypot"
//
//	Color to gray:
//		plambda lena.png "x[0] x[1] x[2] + + 3 /"
//
//	Pick the blue channel of a RGB image:
//		plambda lena.png "x[2]"
//
//	Swap the blue an green channels of a RGB image (6 equivalent ways):
//		plambda lena.png "x[0] x[2] x[1] join3"
//		plambda lena.png "x[0] x[2] x[1] join join"
//		plambda lena.png "x[0] x[1] x[2] rot join3"
//		plambda lena.png "x[0] x[1] x[2] rot join join"
//		plambda lena.png "x split rot join join"
//		plambda lena.png "x split rot join3"
//
//	Merge the two components of a vector field into a single file
//		plambda x.tiff y.tiff "x y join" > xy.tiff
//
//	Set to 0 the green component of a RGB image
//		plambda lena.png "x[0] 0 x[2] join3"
//
//	Naive Canny filter:
//		cat lena.png | gblur 2 | plambda - "x(1,0) 2 * x(1,1) x(1,-1) + + x(-1,0) 2 * x(-1,1) x(-1,-1) + + - >1 x(0,1) 2 * x(1,1) x(-1,1) + + x(0,-1) 2 * x(1,-1) x(-1,-1) + + - >2 <1 <2 hypot <2 <1 atan2 join" | plambda - "x[0] 4 > >1 x[1] fabs pi 4 / > x[1] fabs pi 4 / 3 * < * >2 x[1] fabs pi 4 / < x[1] fabs pi 4 / 3 * > + >3 x[0] x[0](0,1) > x[0] x[0](0,-1) > * >4 x[0] x[0](1,0) > x[0] x[0](-1,0) > * >5 <1 <3 <5 * * <1 <2 <4 * * + x[0] *" | qauto | display
//
//	Anti-Lalpacian (solve Poisson equation):
//		cat lena.png | fft 1 | plambda - "x  :I :I * :J :J * + / -1 *" | fft -1 | qauto | display
//
//
//
// TODO: 2x2 and 3x3 matrix multiplication
// TODO: admit images of different number of channels
// TODO: implement shunting-yard algorithm to admit infix notation
// TODO: handle 3D and nD images
// TODO: merge colonvars and magicvars (the only difficulty lies in naming)


#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PLAMBDA_MAX_TOKENS 2049
#define PLAMBDA_MAX_VARLEN 0x100
#define PLAMBDA_MAX_PIXELDIM 0x100
#define PLAMBDA_MAX_MAGIC 42


#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif
#ifndef FORJ
#define FORJ(n) for(int j=0;j<(n);j++)
#endif
#ifndef FORK
#define FORK(n) for(int k=0;k<(n);k++)
#endif
#ifndef FORL
#define FORL(n) for(int l=0;l<(n);l++)
#endif

//#include "fragments.c"
#include "fail.c"
#include "xmalloc.c"
#include "random.c"
#include "colorcoords.c"
#include "parsenumbers.c"


#define PLAMBDA_CONSTANT 0 // numeric constant
#define PLAMBDA_SCALAR 1   // pixel component
#define PLAMBDA_VECTOR 2   // whole pixel
#define PLAMBDA_OPERATOR 3 // function
#define PLAMBDA_COLONVAR 4 // colon-type variable
#define PLAMBDA_STACKOP 5  // stack operator
#define PLAMBDA_VARDEF 6   // register variable definition (hacky)
#define PLAMBDA_MAGIC 7    // "magic" modifier (requiring cached global data)

// local functions
static double sum_two_doubles      (double a, double b) { return a + b; }
static double substract_two_doubles(double a, double b) { return a - b; }
static double multiply_two_doubles (double a, double b) { return a * b; }
static double divide_two_doubles   (double a, double b) {
	if (!b && !a) return 0;
	return a / b;
}
static double logic_g      (double a, double b) { return a > b; }
static double logic_l      (double a, double b) { return a < b; }
static double logic_e      (double a, double b) { return a == b; }
static double logic_ge     (double a, double b) { return a >= b; }
static double logic_le     (double a, double b) { return a <= b; }
static double logic_ne     (double a, double b) { return a != b; }
static double logic_if (double a, double b, double c) { return a ? b : c; }

static double function_isfinite  (double x) { return isfinite(x); }
static double function_isinf     (double x) { return isinf(x);    }
static double function_isnan     (double x) { return isnan(x);    }
static double function_isnormal  (double x) { return isnormal(x); }
static double function_signbit   (double x) { return signbit(x);  }

static double quantize_255 (double x)
{
	int ix = x;
	if (ix < 0) return 0;
	if (ix > 255) return 255;
	return ix;
}

static double quantize_easy(double x, double a, double b)
{
	return quantize_255(255.0*(x-a)/(b-a));
}

static double range(double x, double a, double b)
{
	return (x-a)/(b-a);
}

static double bound_double(double x, double a, double b)
{
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static void from_cartesian_to_polar(float *y, float *x)
{
	y[0] = hypot(x[0], x[1]);
	y[1] = atan2(x[1], x[0]);
}

static void from_polar_to_cartesian(float *y, float *x)
{
	y[0] = x[0] * cos(x[1]);
	y[1] = x[0] * sin(x[1]);
}

static void complex_product(float *xy, float *x, float *y)
{
	xy[0] = x[0]*y[0] - x[1]*y[1];
	xy[1] = x[0]*y[1] + x[1]*y[0];
}

// table of all functions (local and from math.h)
struct predefined_function {
	void (*f)(void);
	char *name;
	int nargs;
	float value;
} global_table_of_predefined_functions[] = {
#define REGISTER_FUNCTION(x,n) {(void(*)(void))x, #x, n, 0}
#define REGISTER_FUNCTIONN(x,xn,n) {(void(*)(void))x, xn, n, 0}
	REGISTER_FUNCTION(acos,1),
	REGISTER_FUNCTION(acosh,1),
	REGISTER_FUNCTION(asin,1),
	REGISTER_FUNCTION(asinh,1),
	REGISTER_FUNCTION(atan,1),
	REGISTER_FUNCTION(atanh,1),
	REGISTER_FUNCTION(cbrt,1),
	REGISTER_FUNCTION(ceil,1),
	REGISTER_FUNCTION(cos,1),
	REGISTER_FUNCTION(cosh,1),
	REGISTER_FUNCTION(erf,1),
	REGISTER_FUNCTION(erfc,1),
	REGISTER_FUNCTION(exp,1),
	REGISTER_FUNCTION(exp2,1),
	REGISTER_FUNCTION(expm1,1),
	REGISTER_FUNCTION(fabs,1),
	REGISTER_FUNCTION(floor,1),
	REGISTER_FUNCTION(lgamma,1),
	REGISTER_FUNCTION(log,1),
	REGISTER_FUNCTION(log10,1),
	REGISTER_FUNCTION(log1p,1),
	REGISTER_FUNCTION(log2,1),
	REGISTER_FUNCTION(logb,1),
	REGISTER_FUNCTION(nearbyint,1),
	REGISTER_FUNCTION(rint,1),
	REGISTER_FUNCTION(round,1),
	REGISTER_FUNCTION(sin,1),
	REGISTER_FUNCTION(sinh,1),
	REGISTER_FUNCTION(sqrt,1),
	REGISTER_FUNCTION(tan,1),
	REGISTER_FUNCTION(tanh,1),
	REGISTER_FUNCTION(tgamma,1),
	REGISTER_FUNCTION(trunc,1),
	REGISTER_FUNCTION(atan2,2),
	REGISTER_FUNCTION(copysign,2),
	REGISTER_FUNCTION(fdim,2),
	REGISTER_FUNCTION(fma,2),
	REGISTER_FUNCTION(fmax,2),
	REGISTER_FUNCTION(fmin,2),
	REGISTER_FUNCTION(fmod,2),
	REGISTER_FUNCTION(hypot,2),
	REGISTER_FUNCTION(ldexp,2),
	REGISTER_FUNCTION(nextafter,2),
	REGISTER_FUNCTION(nexttoward,2),
	REGISTER_FUNCTION(pow,2),
	REGISTER_FUNCTION(remainder,2),
	REGISTER_FUNCTIONN(quantize_255,"q255",1),
	REGISTER_FUNCTIONN(quantize_easy,"qe",3),
	REGISTER_FUNCTIONN(range,"range",3),
	REGISTER_FUNCTIONN(bound_double,"bound",3),
	REGISTER_FUNCTIONN(pow,"^",2),
	REGISTER_FUNCTIONN(sum_two_doubles,"+",2),
	REGISTER_FUNCTIONN(logic_g,">",2),
	REGISTER_FUNCTIONN(logic_l,"<",2),
	REGISTER_FUNCTIONN(logic_e,"=",2),
	REGISTER_FUNCTIONN(logic_ge,">=",2),
	REGISTER_FUNCTIONN(logic_le,"<=",2),
	REGISTER_FUNCTIONN(logic_ne,"!=",2),
	REGISTER_FUNCTIONN(logic_if,"if",3),
	REGISTER_FUNCTIONN(function_isfinite,"isfinite",1),
	REGISTER_FUNCTIONN(function_isinf,"isinf",1),
	REGISTER_FUNCTIONN(function_isnan,"isnan",1),
	REGISTER_FUNCTIONN(function_isnormal,"isnormal",1),
	REGISTER_FUNCTIONN(function_signbit,"signbit",1),
	REGISTER_FUNCTIONN(divide_two_doubles,"/",2),
	REGISTER_FUNCTIONN(multiply_two_doubles,"*",2),
	REGISTER_FUNCTIONN(substract_two_doubles,"-",2),
	REGISTER_FUNCTIONN(random_uniform,"randu",-1),
	REGISTER_FUNCTIONN(random_normal,"randn",-1),
	REGISTER_FUNCTIONN(random_normal,"randg",-1),
	REGISTER_FUNCTIONN(random_cauchy,"randc",-1),
	REGISTER_FUNCTIONN(random_laplace,"randl",-1),
	REGISTER_FUNCTIONN(random_exponential,"rande",-1),
	REGISTER_FUNCTIONN(random_pareto,"randp",-1),
	REGISTER_FUNCTIONN(random_raw,"rand",-1),
	REGISTER_FUNCTIONN(from_cartesian_to_polar,"topolar", -2),
	REGISTER_FUNCTIONN(from_polar_to_cartesian,"frompolar", -2),
	REGISTER_FUNCTIONN(complex_product,"cprod", -3),
	//REGISTER_FUNCTIONN(rgb2hsv,"rgb2hsv",3),
	//REGISTER_FUNCTIONN(hsv2rgb,"rgb2hsv",3),
#undef REGISTER_FUNCTION
	{NULL, "pi", 0, M_PI},
#ifdef M_E
#define REGISTER_CONSTANT(x) {NULL, #x, 0, x}
	REGISTER_CONSTANT(M_E),
	REGISTER_CONSTANT(M_LOG2E),
	REGISTER_CONSTANT(M_LOG10E),
	REGISTER_CONSTANT(M_LN2),
	REGISTER_CONSTANT(M_LN10),
	REGISTER_CONSTANT(M_PI),
	REGISTER_CONSTANT(M_PI_2),
	REGISTER_CONSTANT(M_PI_4),
	REGISTER_CONSTANT(M_1_PI),
	REGISTER_CONSTANT(M_2_PI),
	REGISTER_CONSTANT(M_2_SQRTPI),
	REGISTER_CONSTANT(M_SQRT2),
	REGISTER_CONSTANT(M_SQRT1_2),
#undef REGISTER_CONSTANT
#endif
};


struct plambda_token {
	int type;
	float value;         // if type==constant, value
	int index;           // if type==variable, its index
	                     // if type==operator, its index
	int component;       // if type==variable, index of selected component
	int displacement[2]; // if type==variable, relative displacement
	int colonvar;        // if type==colon, the letter

	char *tmphack;       // temporary place for storing the unsorted index
};

struct collection_of_varnames {
	int n;
	char *t[PLAMBDA_MAX_TOKENS];
};

struct plambda_program {
	int n;
	struct plambda_token t[PLAMBDA_MAX_TOKENS];
	struct collection_of_varnames var[1];
	int regn[10]; // registers
	float regv[10][PLAMBDA_MAX_PIXELDIM];
};


static float apply_function(struct predefined_function *f, float *v)
{
	switch(f->nargs) {
	case 0: return f->value;
	case 1: return ((double(*)(double))(f->f))(v[0]);
	case 2: return ((double(*)(double,double))f->f)(v[1], v[0]);
	case 3: return ((double(*)(double,double,double))f->f)(v[2],v[1],v[0]);
	case -1: return ((double(*)())(f->f))();
	default: fail("bizarre");
	}
	//return 0;
}

static int symmetrize_index_inside(int i, int m)
{
	assert( i >= 0 && i < m);
	int r = 0;
	if (i > m/2) r = i-m;
	if (i < m/2) r = i;
	return r;
}

// the value of colon variables depends on the position within the image
static float eval_colonvar(int w, int h, int i, int j, int c)
{
	switch(c) {
	case 'i': return i;
	case 'j': return j;
	case 'w': return w;
	case 'h': return h;
	case 'n': return w*h;
	case 'x': return (2.0/(w-1))*i - 1;
	case 'y': return (2.0/(h-1))*j - 1;
	case 'r': return hypot((2.0/(h-1))*j-1,(2.0/(w-1))*i-1);
	case 't': return atan2((2.0/(h-1))*j-1,(2.0/(w-1))*i-1);
	case 'I': return symmetrize_index_inside(i,w);
	case 'J': return symmetrize_index_inside(j,w);
	case 'W': return w/(2*M_PI);
	case 'H': return h/(2*M_PI);
	default: fail("unrecognized colonvar \":%c\"", c);
	}
}

struct image_stats {
	bool init_simple, init_vsimple, init_ordered, init_vordered;
	float scalar_min, scalar_max, scalar_avg, scalar_med;
	//float vector_cmin[PLAMBDA_MAX_PIXELDIM];  // component-wise min
	float vector_n1min[PLAMBDA_MAX_PIXELDIM]; // exemplar with min L1
	float vector_n2min[PLAMBDA_MAX_PIXELDIM]; // exemplar with min L2
	float vector_nimin[PLAMBDA_MAX_PIXELDIM]; // exemplar with min Linf
	//float vector_cmax[PLAMBDA_MAX_PIXELDIM];
	float vector_n1max[PLAMBDA_MAX_PIXELDIM];
	float vector_n2max[PLAMBDA_MAX_PIXELDIM];
	float vector_nimax[PLAMBDA_MAX_PIXELDIM];
	float vector_avg[PLAMBDA_MAX_PIXELDIM];
	float vector_med[PLAMBDA_MAX_PIXELDIM];
	bool init_csimple, init_cordered;
	float component_min[PLAMBDA_MAX_PIXELDIM];
	float component_max[PLAMBDA_MAX_PIXELDIM];
	float component_avg[PLAMBDA_MAX_PIXELDIM];
	float component_med[PLAMBDA_MAX_PIXELDIM];
	float *sorted_samples, *sorted_components[PLAMBDA_MAX_PIXELDIM];
};

struct linear_statistics {
	float min, max, avg, avgnz;
	int n, rns, rnz, nnan, ninf;
};

static void compute_linstats(struct linear_statistics *s,
		float *x, int n, int stride, int offset)
{
	int rns = 0, rnz = 0, nnan = 0, ninf = 0;
	float min = INFINITY, max = -INFINITY;
	long double avg = 0, avgnz = 0;
	for (int i = 0; i < n; i++) {
		float y = x[i*stride + offset];
		if (isnan(y)) {
			nnan += 1;
			continue;
		}
		if (!isfinite(y)) ninf += 1;
		if (y < min) min = y;
		if (y > max) max = y;
		avg += y;
		rns += 1;
		if (y) {
			avgnz += y;
			rnz += 1;
		}
	}
	avg /= rns; avgnz /= rnz;
	s->min=min; s->max=max; s->avg=avg; s->avgnz=avgnz;
	s->n=n; s->rns=rns; s->rnz=rnz; s->nnan=nnan; s->ninf=ninf;
}

static void compute_simple_sample_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_simple) return;
	if (w*h > 1) s->init_simple = true;
	struct linear_statistics ls[1];
	compute_linstats(ls, x, w*h*pd, 1, 0);
	s->scalar_min = ls->min;
	s->scalar_max = ls->max;
	s->scalar_avg = ls->avg;
}

static void compute_simple_component_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_csimple) return;
	if (w*h > 1) s->init_csimple = true;
	for (int l = 0; l < pd; l++)
	{
		struct linear_statistics ls[1];
		compute_linstats(ls, x, w*h, pd, l);
		s->component_min[l] = ls->min;
		s->component_max[l] = ls->max;
		s->component_avg[l] = ls->avg;
	}
}

static float euclidean_norm_of_float_vector(const float *x, int n)
{
	if (n == 1) return fabs(x[0]);
	else {
		float r = 0;
		for (int i = 0; i < n; i++)
			r = hypot(r, x[i]);
		return r;
	}
}

static void compute_simple_vector_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_vsimple) return;
	if (w*h > 1) s->init_vsimple = true;
	int np = w * h, rnp = 0;
	float minpixel = INFINITY, maxpixel = -INFINITY;
	long double avgpixel[PLAMBDA_MAX_PIXELDIM] = {0}, avgnorm = 0;
	int minidx=-1, maxidx=-1;
	for (int j = 0; j < pd; j++)
		avgpixel[j] = 0;
	for (int i = 0; i < np; i++)
	{
		float xnorm = euclidean_norm_of_float_vector(x + pd*i, pd);
		if (isnan(xnorm)) continue;
		if (xnorm < minpixel) { minidx = i; minpixel = xnorm; }
		if (xnorm > maxpixel) { maxidx = i; maxpixel = xnorm; }
		for (int j = 0; j < pd; j++)
			avgpixel[j] += x[pd*i+j];
		avgnorm += xnorm;
		rnp += 1;
	}
	//assert(rnp);
	avgnorm /= rnp;
	long double mipi[pd], mapi[pd];
	for (int j = 0; j < pd; j++) {
		mipi[j] = x[minidx*pd+j];
		mapi[j] = x[maxidx*pd+j];
		avgpixel[j] /= rnp;
	}
	FORI(pd) s->vector_n2min[i] = mipi[i];
	FORI(pd) s->vector_n2max[i] = mapi[i];
	FORI(pd) s->vector_avg[i] = avgpixel[i];
	//setnumber(p, "error", avgnorm);
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

static void compute_ordered_sample_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_ordered) return;
	if (w*h > 1) s->init_ordered = true;
	int ns = w * h * pd;
	s->sorted_samples = xmalloc(ns*sizeof(float));
	FORI(ns) s->sorted_samples[i] = x[i];
	qsort(s->sorted_samples, ns, sizeof(float), compare_floats);
	s->scalar_med = s->sorted_samples[ns/2];
}

static void compute_ordered_component_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_cordered) return;
	if (w*h > 1) s->init_cordered = true;
	int ns = w * h;
	float *t = xmalloc(pd*ns*sizeof(float));
	for (int l = 0; l < pd; l++)
	{
		s->sorted_components[l] = t + l*ns;
		FORI(ns) s->sorted_components[l][i] = x[i*pd+l];
		qsort(s->sorted_components[l],ns,sizeof(float),compare_floats);
		s->component_med[l] = s->sorted_components[l][ns/2];
	}
}

static void compute_ordered_vector_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	fail("ordered vector stats not implemented");
	(void)x;
	(void)w;
	(void)h;
	(void)pd;
	// there is some bizarre trickery waiting to be coded in here
	s->init_vordered = true;
}

static int bound(int a, int x, int b)
{
	if (b < a) return bound(b, x, a);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}


// the value of magic variables depends on some globally cached data
static int eval_magicvar(float *out, int magic, int img_index, int comp, int qq,
		float *x, int w, int h, int pd) // only needed on the first run
{
	// XXX WARNING : global variables here (leading to non-re-entrant code)
	static bool initt = false;
	//static struct image_stats *t = 0;
	static struct image_stats t[PLAMBDA_MAX_MAGIC];
	if (!initt) {
		//t = xmalloc(PLAMBDA_MAX_MAGIC * sizeof*t);
		for (int i = 0; i < PLAMBDA_MAX_MAGIC; i++) {
			t[i].init_simple = false;
			t[i].init_ordered = false;
			t[i].init_vsimple = false;
			t[i].init_vordered = false;
			t[i].init_csimple = false;
			t[i].init_cordered = false;
		}
		initt = true;
	}
	//fprintf(stderr, "magic=%c index=%d comp=%d\n",magic,img_index,comp);

	if (img_index >= PLAMBDA_MAX_MAGIC)
		fail("%d magic images is too much for me!", PLAMBDA_MAX_MAGIC);

	struct image_stats *ti = t + img_index;

	if (magic=='i' || magic=='a' || magic=='v') {
		if (comp < 0) { // use all samples
			compute_simple_sample_stats(ti, x, w, h, pd);
			switch(magic) {
				case 'i': *out = ti->scalar_min; break;
				case 'a': *out = ti->scalar_max; break;
				case 'v': *out = ti->scalar_avg; break;
				default: fail("this can not happen");
			}
			return 1;
		} else { // use samples from the specified component
			compute_simple_component_stats(ti, x, w, h, pd);
			switch(magic) {
				case 'i': *out = ti->component_min[comp]; break;
				case 'a': *out = ti->component_max[comp]; break;
				case 'v': *out = ti->component_avg[comp]; break;
				default: fail("this can not happen");
			}
			return 1;
		}
	} else if (magic=='I' || magic=='A' || magic=='V') {
		compute_simple_vector_stats(ti, x, w, h, pd);
		switch(magic) {
		case 'I': FORI(pd) out[i] = ti->vector_n2min[i]; break;
		case 'A': FORI(pd) out[i] = ti->vector_n2max[i]; break;
		case 'V': FORI(pd) out[i] = ti->vector_avg[i]; break;
		default: fail("this can not happen");
		}
		return pd;
	} else if (magic=='Y' || magic=='E') {
		compute_simple_component_stats(ti, x, w, h, pd);
		switch(magic) {
		case 'Y': FORI(pd) out[i] = ti->component_min[i]; break;
		case 'E': FORI(pd) out[i] = ti->component_max[i]; break;
		default: fail("this can not happen");
		}
		return pd;
	} else if (magic == 'm' || magic == 'q') {
		if (comp < 0) { // use all samples
			compute_ordered_sample_stats(ti, x, w, h, pd);
			if (magic == 'm') {
				*out = ti->scalar_med;
				return 1;
			}
			if (magic == 'q') {
				int qpos = round(qq*w*h*pd/100.0);
				qpos = bound(0, qpos, w*h*pd-1);
				*out = ti->sorted_samples[qpos];
				return 1;
			}
		} else {
			compute_ordered_component_stats(ti, x, w, h, pd);
			if (magic == 'm') {
				*out = ti->component_med[comp];
				return 1;
			}
			if (magic == 'q') {
				int qpos = round(qq*w*h/100.0);
				qpos = bound(0, qpos, w*h-1);
				*out = ti->sorted_components[comp][qpos];
				return 1;
			}
		}
	} else
		fail("magic of kind '%c' is not yed implemented", magic);

	return 0;
}


// if the token resolves to a numeric constant, store it in *x and return true
// otherwise, return false
// if trailing characters are ignored, print a warning message
static bool token_is_number(float *x, const char *t)
{
	char *endptr;
	*x = strtof(t, &endptr);
	if (endptr == t) return false;
	if (*endptr != '\0')
		fprintf(stderr, "TOKEN "
				"WARNING: trailing characters (\"%s\") "
				"ignored "
			       	"in numeric constant\n", endptr);
	return true;
}

// if token is colonvar, return the id
// otherwise, return zero
static int token_is_colonvar(const char *t)
{
	if (t[0] != ':') return 0;
	if (isalpha(t[1]) && t[2]=='\0') return t[1];
	return 0;
}

// if token is a variable definition, return the index
// otherwise, return zero
static int token_is_vardef(const char *t)
{
	if (t[0]=='>' && isdigit(t[1]) && t[1]>'0' && t[2]=='\0')
		return t[1] - '0';
	if (t[0]=='<' && isdigit(t[1]) && t[1]>'0' && t[2]=='\0')
		return -(t[1] - '0');
	return 0;
}


#define PLAMBDA_STACKOP_NO 0
#define PLAMBDA_STACKOP_DEL 1
#define PLAMBDA_STACKOP_DUP 2
#define PLAMBDA_STACKOP_VSPLIT 3
#define PLAMBDA_STACKOP_VMERGE 4
#define PLAMBDA_STACKOP_ROT 5
#define PLAMBDA_STACKOP_VMERGE3 6
#define PLAMBDA_STACKOP_VMERGEALL 7
#define PLAMBDA_STACKOP_HSV2RGB 8
#define PLAMBDA_STACKOP_RGB2HSV 9
#define PLAMBDA_STACKOP_NMERGE 10
#define PLAMBDA_STACKOP_INTERLEAVE 11
#define PLAMBDA_STACKOP_DEINTERLEAVE 12

// if token is a stack operation, return its id
// otherwise, return zero
static int token_is_stackop(const char *t)
{
	if (0 == strcmp(t, "del")) return PLAMBDA_STACKOP_DEL;
	if (0 == strcmp(t, "dup")) return PLAMBDA_STACKOP_DUP;
	if (0 == strcmp(t, "rot")) return PLAMBDA_STACKOP_ROT;
	if (0 == strcmp(t, "split")) return PLAMBDA_STACKOP_VSPLIT;
	if (0 == strcmp(t, "merge")) return PLAMBDA_STACKOP_VMERGE;
	if (0 == strcmp(t, "join")) return PLAMBDA_STACKOP_VMERGE;
	if (0 == strcmp(t, "merge3")) return PLAMBDA_STACKOP_VMERGE3;
	if (0 == strcmp(t, "join3")) return PLAMBDA_STACKOP_VMERGE3;
	if (0 == strcmp(t, "mergeall")) return PLAMBDA_STACKOP_VMERGEALL;
	if (0 == strcmp(t, "joinall")) return PLAMBDA_STACKOP_VMERGEALL;
	if (0 == strcmp(t, "hsv2rgb")) return PLAMBDA_STACKOP_HSV2RGB;
	if (0 == strcmp(t, "rgb2hsv")) return PLAMBDA_STACKOP_RGB2HSV;
	if (0 == strcmp(t, "njoin")) return PLAMBDA_STACKOP_NMERGE;
	if (0 == strcmp(t, "nmerge")) return PLAMBDA_STACKOP_NMERGE;
	if (0 == strcmp(t, "interleave")) return PLAMBDA_STACKOP_INTERLEAVE;
	if (0 == strcmp(t, "deinterleave")) return PLAMBDA_STACKOP_DEINTERLEAVE;
	return 0;
}

// if the token is a valid word, return its length
//         and if the token is followed by modifiers, fill *endptr
// otherwise, return zero
static int token_is_word(const char *t, const char **endptr)
{
	*endptr = NULL;
	if ((*t=='+'||*t=='-'||*t=='/'||*t=='^'||*t=='*'||*t=='>'||*t=='<'||*t=='=')&&t[1]=='\0')
		return 1;
	if (!isalpha(t[0])) {
		return 0;
	}
	int n = 1;
	while (t[n]) {
		if  (!isalnum(t[n])) {
			*endptr = t+n;
			return (t[n]=='(' || t[n]=='[' || t[n]=='%') ? n : 0;
		}
		n += 1;
	}
	return n;
}

static int word_is_predefined(const char *id)
{
	int n = sizeof(global_table_of_predefined_functions)/
		sizeof(global_table_of_predefined_functions[0]);
	struct predefined_function *r = global_table_of_predefined_functions;
	FORI(n)
		if (0 == strcmp(r[i].name, id))
			return i;
	return -1;
}

// fills the modifiers with their defined values, otherwise with the default
static void parse_modifiers(const char *mods,
		int *ocomp, int *odx, int *ody, int *omagic)
{
	*ocomp = -1;
	*odx = 0;
	*ody = 0;
	*omagic = 0;
	int comp, dx, dy;
	char magic;
	// NOTE: the order of the following conditions is important for a
	// correct parsing
	if (!mods) {
		return;
	} else if (3 == sscanf(mods, "[%d](%d,%d)", &comp, &dx, &dy)) {
		*odx = dx;
		*ody = dy;
		*ocomp = comp;
	 	return;
	} else if (3 == sscanf(mods, "(%d,%d)[%d]", &dx, &dy, &comp)) {
		*odx = dx;
		*ody = dy;
		*ocomp = comp;
	 	return;
	} else if (2 == sscanf(mods, "(%d,%d)", &dx, &dy)) {
		*odx = dx;
		*ody = dy;
	 	return;
	} else if (2 == sscanf(mods, "%%%c%d", &magic, &dx)) {
		*omagic = magic;
		*odx = dx;
		return;
	} else if (1 == sscanf(mods, "%%%c", &magic)) {
		*omagic = magic;
		return;
	} else if (3 == sscanf(mods, "[%d]%%%c%d", &comp, &magic, &dx)) {
		*omagic = magic;
		*ocomp = comp;
		*odx = dx;
		return;
	} else if (2 == sscanf(mods, "[%d]%%%c", &comp, &magic)) {
		*omagic = magic;
		*ocomp = comp;
		return;
	} else if (1 == sscanf(mods, "[%d]", &comp)) {
		*ocomp = comp;
		return;
	}
}

static void collection_of_varnames_init(struct collection_of_varnames *x)
{
	x->n = 0;
}

static int collection_of_varnames_find(struct collection_of_varnames *x,
		const char *s)
{
	FORI(x->n)
		if (0 == strcmp(s, x->t[i]))
			return i;
	return -1;
}

static char *collection_of_varnames_add(struct collection_of_varnames *x,
		const char *s)
{
	char *r;
	int i = collection_of_varnames_find(x, s);
	if (i < 0) {
		if (x->n+1 >= PLAMBDA_MAX_TOKENS)
			fail("caca");
		r = xmalloc(1+strlen(s));
		strcpy(r, s);
		x->t[x->n] = r;
		x->n += 1;
	} else {
		r = x->t[i];
	}
	return r;
}

static void collection_of_varnames_end(struct collection_of_varnames *x)
{
	FORI(x->n)
		free(x->t[i]);
	x->n = 0;
}


static int strcmp_for_qsort(const void *aa, const void *bb)
{
	const char **a = (const char **)aa;
	const char **b = (const char **)bb;
	return strcmp(*a, *b);
}

static void collection_of_varnames_sort(struct collection_of_varnames *x)
{
	qsort(x->t, x->n, sizeof*x->t, strcmp_for_qsort);
}

// this function takes a string which contains one token,
// and compiles the corresponding info into p->t[p->n]
//
// TODO (maybe): split token identification from info gathering
// (this will produce longer code but shorter functions)
static void process_token(struct plambda_program *p, const char *tokke)
{
	char tok[1+strlen(tokke)];             // the string of the token
	strcpy(tok, tokke);
	struct plambda_token *t = p->t + p->n; // the compiled token

	int tok_id;
	const char *tok_end;

	float x;
	if (token_is_number(&x, tok)) {
		t->type = PLAMBDA_CONSTANT;
		t->value = x;
		goto endtok;
	}

	if ((tok_id = token_is_colonvar(tok))) {
		t->type = PLAMBDA_COLONVAR;
		t->colonvar = tok_id;
		goto endtok;
	}

	if ((tok_id = token_is_stackop(tok))) {
		t->type = PLAMBDA_STACKOP;
		t->index = tok_id;
		goto endtok;
	}

	if ((tok_id = token_is_vardef(tok))) {
		t->type = PLAMBDA_VARDEF;
		t->index = tok_id;
		goto endtok;
	}

	if ((token_is_word(tok, &tok_end)))
	{
		int idx = word_is_predefined(tok);
		if (idx < 0) {
			char varname[PLAMBDA_MAX_VARLEN+1];
			int varlen = strlen(tok);
			if (tok_end) varlen = tok_end-tok;
			if (varlen >= PLAMBDA_MAX_VARLEN)
				varlen = PLAMBDA_MAX_VARLEN;
			FORI(varlen) varname[i] = tok[i];
			varname[varlen] = '\0';
			int comp, disp[2], magic;
			t->tmphack =collection_of_varnames_add(p->var, varname);
			parse_modifiers(tok_end, &comp, disp, disp+1, &magic);
			t->type = comp<0 ? PLAMBDA_VECTOR : PLAMBDA_SCALAR;
			if (magic) {
				t->type = PLAMBDA_MAGIC;
				t->colonvar = magic;
			}
			t->component = comp;
			t->displacement[0] = disp[0];
			t->displacement[1] = disp[1];
		} else {
			//struct predefined_function *f =
			//	global_table_of_predefined_functions + idx;
			t->type = PLAMBDA_OPERATOR;
			t->index = idx;
		}
		goto endtok;
	}

endtok:
	p->n += 1;
}

// this function updates the indexes of a
// collection of variables which is sorted in alphabetical order
static void unhack_varnames(struct plambda_program *p)
{
	FORI(p->n)
	{
		struct plambda_token *t = p->t + i;
		if (t->type == PLAMBDA_SCALAR || t->type == PLAMBDA_VECTOR
				|| t->type == PLAMBDA_MAGIC)
		{
			t->index = collection_of_varnames_find(p->var,
								t->tmphack);
			if (t->index < 0)
				fail("unexpected bad variable \"%s\"",
								t->tmphack);
		}
	}
}

static void plambda_compile_program(struct plambda_program *p, const char *str)
{
	char s[1+strlen(str)];
	strcpy(s, str);
	char *spacing = " \n\t";

	FORI(10) p->regn[i] = 0;

	collection_of_varnames_init(p->var);
	p->n = 0;
	int n = 0;
	char *tok = strtok(s, spacing);
	while (tok) {
		//fprintf(stderr, "token[%d] = %s\n", n, tok);
		process_token(p, tok);
		tok = strtok(NULL, spacing);
		n += 1;
	}

	collection_of_varnames_sort(p->var);

	// the "sort" above does not update the variable indices
	// the following function updates them
	unhack_varnames(p);
}

static const char *arity(struct predefined_function *f)
{
	switch(f->nargs) {
	case 0: return "0-ary";
	case 1: return "unary";
	case 2: return "binary";
	case 3: return "ternary";
	case -1: return "strange";
	case -2: return "strange2";
	case -3: return "strange3";
	default: return "unrecognized";
	}
}

static void print_compiled_program(struct plambda_program *p)
{
	fprintf(stderr, "COMPILED PROGRAM OF %d TOKENS:\n", p->n);
	FORI(p->n) {
		struct plambda_token *t = p->t + i;
		fprintf(stderr, "TOKEN[%d]: ", i);
		if (t->type == PLAMBDA_CONSTANT)
			fprintf(stderr, "constant %g", t->value);
		if (t->type == PLAMBDA_COLONVAR)
			fprintf(stderr, "colonvar \"%c\"", t->colonvar);
		if (t->type == PLAMBDA_VECTOR) {
			fprintf(stderr, "variable vector %d \"%s\"",
					t->index, p->var->t[t->index]);
			fprintf(stderr, ", displacement (%d,%d)",
					t->displacement[0], t->displacement[1]);
		}
		if (t->type == PLAMBDA_SCALAR) {
			fprintf(stderr, "variable scalar %d \"%s\"",
					t->index, p->var->t[t->index]);
			fprintf(stderr, ", displacement (%d,%d)",
					t->displacement[0], t->displacement[1]);
			fprintf(stderr, ", component %d", t->component);
		}
		if (t->type == PLAMBDA_OPERATOR) {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			fprintf(stderr, "%s operator %s", arity(f), f->name);
		}
		if (t->type == PLAMBDA_STACKOP)
			fprintf(stderr, "stack manipulation");
		if (t->type == PLAMBDA_VARDEF) {
			fprintf(stderr, "register variable %s %d",
					t->index<0 ? "read" : "write",
					abs(t->index));
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "The program uses %d variables:\n", p->var->n);
	FORI(p->var->n)
		fprintf(stderr, "VARIABLE[%d] = \"%s\"\n", i, p->var->t[i]);
}


// stack of vectorial (or possibly scalar) values
struct value_vstack {
	int n;
	int d[PLAMBDA_MAX_TOKENS];
	float t[PLAMBDA_MAX_TOKENS][PLAMBDA_MAX_PIXELDIM];
};

static int vstack_pop_vector(float *val, struct value_vstack *s)
{
	if (s->n > 0) {
		s->n -= 1;
		int d = s->d[s->n];
		if (val) FORI(d) val[i] = s->t[s->n][i];
		return d;
	} else fail("popping from empty stack");
}

static void vstack_push_vector(struct value_vstack *s, float *v, int n)
{
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		s->d[s->n] = n;
		FORI(n)
			s->t[s->n][i] = v[i];
		s->n += 1;
	} else fail("full stack");
}

static void vstack_push_scalar(struct value_vstack *s, float x)
{
	vstack_push_vector(s, &x, 1);
	//if (s->n+1 < PLAMBDA_MAX_TOKENS) {
	//	s->d[s->n] = 1;
	//	s->t[s->n][0] = x;
	//	s->n += 1;
	//} else fail("full stack");
}

//static void vstack_print(FILE *f, struct value_vstack *s)
//{
//	FORI(s->n) {
//		fprintf(f, "STACK[%d/%d]: {%d}", 1+i, s->n, s->d[i]);
//		FORJ(s->d[i])
//			fprintf(f, " %g", s->t[i][j]);
//		fprintf(f, "\n");
//	}
//}



// XXX TODO: refactor the following strange cases into the general setup

// a function that takes no arguments but must be called nonetheless
static void treat_strange_case(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -1);
	float r = apply_function(f, NULL);
	vstack_push_vector(s, &r, 1);
}

// pop a 2-vector x from the stack and push f(x) as a 2-vector
static void treat_strange_case2(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -2);
	float v[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int n = vstack_pop_vector(v, s);
	if (n != 2) fail("function \"%s\" requires a 2-vector", f->name);
	((void(*)(float*,float*))(f->f))(r, v);
	vstack_push_vector(s, r, n);
}

#ifndef ODDP
#define ODDP(x) ((x)&1)
#endif
#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

// pop two 2-vectors x,y from the stack and push f(x,y) as a 2-vector
// NOTE: works also on vectors of even length, doing the "obvious" thing
static void treat_strange_case3(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -3);
	float a[PLAMBDA_MAX_PIXELDIM];
	float b[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int na = vstack_pop_vector(a, s);
	int nb = vstack_pop_vector(b, s);
	if (ODDP(na)) fail("function \"%s\" requires 2-vectors", f->name);
	if (ODDP(nb)) fail("function \"%s\" requires 2-vectors", f->name);
	void (*ff)(float*,float*,float*) = (void(*)(float*,float*,float*))f->f;
	int ca = na / 2;
	int cb = nb / 2;
	int n = 0;
	if (ca == cb) {
		n = na;
		for (int i = 0; i < ca; i++)
			ff(r+2*i, a+2*i, b+2*i);
	} else if (ca == 1) {
		n = nb;
		for (int i = 0; i < cb; i++)
			ff(r+2*i, a, b+2*i);
	} else if (cb == 1) {
		n = na;
		for (int i = 0; i < ca; i++)
			ff(r+2*i, a+2*i, b);
	} else fail("function \"%s\" can not operate on lengths (%d,%d)",
			f->name, na, nb);
	assert(n);
	//if (na != nb) fail("this can not happen");
	//((void(*)(float*,float*,float*))(f->f))(r, a, b);
	vstack_push_vector(s, r, n);
}

// pop a 3-vector x from the stack and push f(x) as a 3-vector
static void treat_strange_case4(struct value_vstack *s,
		struct predefined_function *f)
{
	fail("color space conversions not implemented");
}

// this function is complicated because it contains the scalar+vector
// semantics, which is complicated
static void vstack_apply_function(struct value_vstack *s,
					struct predefined_function *f)
{
	if (f->nargs == -1) {treat_strange_case(s,f); return;}
	if (f->nargs == -2) {treat_strange_case2(s,f); return;}
	if (f->nargs == -3) {treat_strange_case3(s,f); return;}
	int d[f->nargs], rd = 1;
	float v[f->nargs][PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	FORI(f->nargs)
		d[i] = vstack_pop_vector(v[i], s);
	// the d[i] which are larger than one must be equal
	FORI(f->nargs)
		if (d[i] > 1) {
			if (rd > 1 && d[i] != rd)
				fail("can not vectorize (%d %d)", rd, d[i]);
			else
				rd = d[i];
		}
	if (rd > 1)
		FORI(f->nargs)
			if (d[i] == 1)
				FORL(rd)
					v[i][l] = v[i][0];
	FORL(rd) {
		float a[f->nargs];
		FORI(f->nargs)
			a[i] = v[i][l];
		r[l] = apply_function(f, a);
	}
	vstack_push_vector(s, r, rd);
}

static void vstack_process_op(struct value_vstack *s, int opid)
{
	switch(opid) {
	case PLAMBDA_STACKOP_DEL:
		vstack_pop_vector(NULL, s);
		break;
	case PLAMBDA_STACKOP_DUP: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		vstack_push_vector(s, x, n);
		vstack_push_vector(s, x, n);
				  }
		break;
	case PLAMBDA_STACKOP_VSPLIT: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		FORI(n)
			vstack_push_scalar(s, x[i]);
				     }
		break;
	case PLAMBDA_STACKOP_VMERGE: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int m = vstack_pop_vector(y, s);
		int n = vstack_pop_vector(x, s);
		if (n+m >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		FORI(m)
			x[n+i] = y[i];
		vstack_push_vector(s, x, n+m);
				     }
		break;
	case PLAMBDA_STACKOP_VMERGE3: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		float z[PLAMBDA_MAX_PIXELDIM];
		int nz = vstack_pop_vector(z, s);
		int ny = vstack_pop_vector(y, s);
		int nx = vstack_pop_vector(x, s);
		if (nx+ny+nz >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		FORI(ny) x[nx+i] = y[i];
		FORI(nz) x[nx+ny+i] = z[i];
		vstack_push_vector(s, x, nx+ny+nz);
				     }
		break;
	case PLAMBDA_STACKOP_HSV2RGB: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (n != 3) fail("hsv2rgb needs a 3-vector");
		double dx[3] = {x[0], x[1], x[2]};
		double dy[3];
		hsv_to_rgb_doubles(dy, dx);
		FORI(3) x[i] = dy[i];
		vstack_push_vector(s, x, 3);
				      }
		break;
	case PLAMBDA_STACKOP_RGB2HSV: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (n != 3) fail("rgb2hsv needs a 3-vector");
		double dx[3] = {x[0], x[1], x[2]};
		double dy[3];
		rgb_to_hsv_doubles(dy, dx);
		FORI(3) x[i] = dy[i];
		vstack_push_vector(s, x, 3);
				      }
		break;
	case PLAMBDA_STACKOP_ROT: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		int m = vstack_pop_vector(y, s);
		vstack_push_vector(s, x, n);
		vstack_push_vector(s, y, m);
		break;
				  }
	//case PLAMBDA_STACKOP_MUL22: {
	//	float a[6];
	//	FORI(6) {
	//		float tmp[PLAMBDA_MAX_PIXELDIM];
	//		if (1 == vstack_pop_vector)
	//	vstack_pop_vector
	//			    }
	case PLAMBDA_STACKOP_VMERGEALL:
		fail("mergeall not implemented");
		break;
	case PLAMBDA_STACKOP_NMERGE: {
		float nn[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(nn, s);
		if (n != 1 || nn[0] < 1 || round(nn[0]) != nn[0]
				|| nn[0] >= PLAMBDA_MAX_PIXELDIM)
			fail("can not nmerge \"%d\" things", nn[0]);
		n = nn[0];
		float x[n][PLAMBDA_MAX_PIXELDIM], y[PLAMBDA_MAX_PIXELDIM];
		int d[n], sdi = 0;
		FORI(n) {
			d[i] = vstack_pop_vector(x[i], s);
			sdi + d[i];
		}
		if (sdi >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		int cx = 0;
		FORI(n) FORJ(d[i]) y[cx++] = x[i][j];
		assert(cx == sdi);
		vstack_push_vector(s, y, sdi);
		break;
				     }
	case PLAMBDA_STACKOP_INTERLEAVE: { // strictly speaking, not a stackop
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (ODDP(n)) fail("can not interleave an odd number %d!", n);
		FORI(n/2) {
			y[2*i] = x[i];
			y[2*i+1] = x[i+n/2];
		}
		vstack_push_vector(s, y, n);
		break;
					 }
	case PLAMBDA_STACKOP_DEINTERLEAVE: { // strictly speaking, not a stackop
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (ODDP(n)) fail("can not deinterleave an odd number %d!", n);
		FORI(n/2) {
			y[i] = x[2*i];
			y[i+n/2] = x[2*i+1];
		}
		vstack_push_vector(s, y, n);
		break;
					 }
	default:
		fail("impossible condition (stackop %d)", opid);
	}
}


#include "getpixel.c"

// returns the dimension of the output
static int run_program_vectorially_at(float *out, struct plambda_program *p,
		float **val, int w, int h, int *pd, int ai, int aj)
{
	struct value_vstack s[1];
	s->n = 0;
	FORI(p->n) {
		struct plambda_token *t = p->t + i;
		switch(t->type) {
		case PLAMBDA_STACKOP:
			vstack_process_op(s, t->index);
			break;
		case PLAMBDA_CONSTANT:
			vstack_push_scalar(s, t->value);
			break;
		case PLAMBDA_COLONVAR: {
			float x = eval_colonvar(w, h, ai, aj, t->colonvar);
			vstack_push_scalar(s, x);
			break;
				       }
		case PLAMBDA_SCALAR: {
			float *img = val[t->index];
			int dai = ai + t->displacement[0];
			int daj = aj + t->displacement[1];
			int cmp = t->component;
			int pdv = pd[t->index];
			float x = getsample_1(img, w, h, pdv, dai, daj, cmp);
			vstack_push_scalar(s, x);
			break;
				     }
		case PLAMBDA_VECTOR: {
			float *img = val[t->index];
			int dai = ai + t->displacement[0];
			int daj = aj + t->displacement[1];
			int pdv = pd[t->index];
			float x[pdv];
			FORL(pdv)
				x[l] = getsample_1(img, w, h, pdv, dai, daj, l);
			vstack_push_vector(s, x, pdv);
				     }
			break;
		case PLAMBDA_OPERATOR: {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			vstack_apply_function(s, f);
				       }
			break;
		case PLAMBDA_VARDEF: {
			int n = abs(t->index);
			if (t->index > 0)
				p->regn[n] = vstack_pop_vector(p->regv[n], s);
			if (t->index < 0)
				vstack_push_vector(s, p->regv[n], p->regn[n]);
				     }
			break;
		case PLAMBDA_MAGIC: {
			int pdv = pd[t->index];
			float *img = val[t->index], x[pdv];
			int rm = eval_magicvar(x, t->colonvar, t->index,
					t->component, t->displacement[0],
					img, w, h, pdv);
			vstack_push_vector(s, x, rm);
				    }
			break;
		default:
			fail("unknown tag type %d", t->type);
		}
	}
	return vstack_pop_vector(out, s);
}

static int eval_dim(struct plambda_program *p, float **val, int *pd)
{
	float result[PLAMBDA_MAX_PIXELDIM];
	int r = run_program_vectorially_at(result, p, val, 1, 1, pd, 0, 0);
	return r;
}

// returns the dimension of the output
static int run_program_vectorially(float *out, int pdmax,
		struct plambda_program *p,
		float **val, int w, int h, int *pd)
{
	int r = 0;
	FORJ(h) FORI(w) {
		float result[pdmax];
		r = run_program_vectorially_at(result, p,val, w,h,pd, i,j);
		assert(r == pdmax);
		FORL(r) {
			setsample_0(out, w, h, pdmax, i, j, l, result[l]);
		}
	}
	return r;
}

//static void shrink_components(float *y, float *x, int n, int ypd, int xpd)
//{
//	assert(ypd <= xpd);
//	FORI(n)
//		FORL(ypd)
//			y[ypd*i + l] = x[xpd*i + l];
//}

#include "smapa.h"
SMART_PARAMETER_SILENT(SRAND,0)
SMART_PARAMETER_SILENT(PLAMBDA_CALC,0)

int main_calc(int c, char *v[])
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s v1 v2 ... \"plambda\"\n", *v);
		//                          0 1  2        c-1
		return EXIT_FAILURE;
	}

	struct plambda_program p[1];
	plambda_compile_program(p, v[c-1]);
	//print_compiled_program(p);

	int n = c - 2, pd[n], pdmax = PLAMBDA_MAX_PIXELDIM;
	if (n != p->var->n)
		fail("the program expects %d variables but %d vectors "
					"were given", p->var->n, n);

	float *x[n];
	FORI(n) x[i] = alloc_parse_floats(pdmax, v[i+1], pd+i);

	FORI(n)
		fprintf(stderr, "calculator correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	float out[pdmax];
	int od = run_program_vectorially_at(out, p, x, 1, 1, pd, 0, 0);

	for (int i = 0; i < od; i++)
		printf("%g%c", out[i], i==(od-1)?'\n':' ');

	return EXIT_SUCCESS;
}

#include "iio.h"
int main_images(int c, char *v[])
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s in1 in2 ... \"plambda\"\n", *v);
		//                          0 1   2         c-1
		return EXIT_FAILURE;
	}

	struct plambda_program p[1];

	plambda_compile_program(p, v[c-1]);
	//print_compiled_program(p);

	int n = c - 2;
	if (n != p->var->n && !(n == 1 && p->var->n == 0))
		fail("the program expects %d variables but %d images "
					"were given", p->var->n, n);
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i+1], w + i, h + i, pd + i);
	FORI(n-1)
		if (w[0] != w[i+1] || h[0] != h[i+1])// || pd[0] != pd[i+1])
			fail("input images size mismatch");

	if (n>1) FORI(n)
		fprintf(stderr, "plambda correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	srand(SRAND());

	int pdreal = eval_dim(p, x, pd);

	float *out = xmalloc(*w * *h * pdreal * sizeof*out);
	int opd = run_program_vectorially(out, pdreal, p, x, *w, *h, pd);
	assert(opd == pdreal);

	iio_save_image_float_vec("-", out, *w, *h, opd);

	FORI(n) free(x[i]);
	free(out);
	collection_of_varnames_end(p->var);

	return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
	int (*f)(int c, char *v[]);
       	f = (PLAMBDA_CALC()>0 || **v=='c' || c==2) ? main_calc : main_images;
	return f(c,v);
}
