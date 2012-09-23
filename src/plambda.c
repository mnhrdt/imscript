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
//
//	Merge the two components of a vector field into a single file
//		plambda x.tiff y.tiff "x y join" > xy.tiff
//
//	Set to 0 the green component of a RGB image
//		plambda lena.png "x[0] 0 x[2] join3"
//
//	Naive Canny filter:
//		cat lena.png | gblur 2 | plambda - "x(1,0) 2 * x(1,1) x(1,-1) + + x(-1,0) 2 * x(-1,1) x(-1,-1) + + - >1 x(0,1) 2 * x(1,1) x(-1,1) + + x(0,-1) 2 * x(1,-1) x(-1,-1) + + - >2 <1 <2 hypot <2 <1 atan2 join" | plambda - "x[0] 4 > >1 x[1] fabs pi 4 / > x[1] fabs pi 4 / 3 * < * >2 x[1] fabs pi 4 / < x[1] fabs pi 4 / 3 * > + >3 x[0] x[0](0,1) > x[0] x[0](0,-1) > * >4 x[0] x[0](1,0) > x[0] x[0](-1,0) > * >5 <1 <3 <5 * * <1 <2 <4 * * + x[0] *" | qauto | d
//
//
//
// TODO: 2x2 and 3x3 matrix multiplication
// TODO: admit images of different number of channels
// TODO: implement shunting-yard algorithm to admit infix notation
// TODO: handle 3D and nD images


#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#define PLAMBDA_MAX_TOKENS 2049
#define PLAMBDA_MAX_VARLEN 0x100
#define PLAMBDA_MAX_PIXELDIM 0x100



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


#define PLAMBDA_CONSTANT 0   // numeric constant
#define PLAMBDA_SCALAR 1     // pixel component
#define PLAMBDA_VECTOR 2     // whole pixel
#define PLAMBDA_OPERATOR 3   // function
#define PLAMBDA_COLONVAR 4   // colon-type variable
#define PLAMBDA_STACKOP 5    // stack operator
#define PLAMBDA_VARDEF 6     // register variable definition (hacky)

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
	REGISTER_FUNCTIONN(pow,"^",2),
	REGISTER_FUNCTIONN(sum_two_doubles,"+",2),
	REGISTER_FUNCTIONN(logic_g,">",2),
	REGISTER_FUNCTIONN(logic_l,"<",2),
	REGISTER_FUNCTIONN(logic_e,"=",2),
	REGISTER_FUNCTIONN(logic_ge,">=",2),
	REGISTER_FUNCTIONN(logic_le,"<=",2),
	REGISTER_FUNCTIONN(logic_ne,"!=",2),
	REGISTER_FUNCTIONN(logic_if,"if",3),
	REGISTER_FUNCTIONN(divide_two_doubles,"/",2),
	REGISTER_FUNCTIONN(multiply_two_doubles,"*",2),
	REGISTER_FUNCTIONN(substract_two_doubles,"-",2),
	REGISTER_FUNCTIONN(random_uniform,"randu",-1),
	REGISTER_FUNCTIONN(random_normal,"randn",-1),
	REGISTER_FUNCTIONN(random_cauchy,"randc",-1),
	REGISTER_FUNCTIONN(random_laplace,"randl",-1),
	REGISTER_FUNCTIONN(random_exponential,"rande",-1),
	REGISTER_FUNCTIONN(random_pareto,"randp",-1),
	//REGISTER_FUNCTIONN(rgb2hsv,"rgb2hsv",3),
	//REGISTER_FUNCTIONN(hsv2rgb,"rgb2hsv",3),
#undef REGISTER_FUNCTION
	{NULL, "pi", 0, M_PI}
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
	default: fail("unrecognized colonvar \":%c\"", c);
	}
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
			return (t[n]=='(' || t[n]=='[') ? n : 0;
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
static void parse_modifiers(const char *mods, int *ocomp, int *odx, int *ody)
{
	*ocomp = -1;
	*odx = 0;
	*ody = 0;
	int comp, dx, dy;
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
			int comp, disp[2];
			t->tmphack =collection_of_varnames_add(p->var, varname);
			parse_modifiers(tok_end, &comp, disp, disp+1);
			t->type = comp<0 ? PLAMBDA_VECTOR : PLAMBDA_SCALAR;
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
		if (t->type == PLAMBDA_SCALAR || t->type == PLAMBDA_VECTOR)
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
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		s->d[s->n] = 1;
		s->t[s->n][0] = x;
		s->n += 1;
	} else fail("full stack");
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

static void treat_strange_case(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -1);
	float r = apply_function(f, NULL);
	vstack_push_vector(s, &r, 1);
}

// this function is complicated because it contains the scalar+vector
// semantics, which is complicated
static void vstack_apply_function(struct value_vstack *s,
					struct predefined_function *f)
{
	if (f->nargs == -1) {treat_strange_case(s,f); return;}
	int d[f->nargs], rd = 1;
	float v[f->nargs][PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	FORI(f->nargs)
		d[i] = vstack_pop_vector(v[i], s);
	FORI(f->nargs)
		// TODO: solve commutativity issue here
		if (d[i] != rd) {
			if (rd > 1)
				fail("can not mix %d- and %d-vectors [%s]\n",
					rd, d[i], f->name);
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
	default:
		fail("impossible condition");
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

static void shrink_components(float *y, float *x, int n, int ypd, int xpd)
{
	assert(ypd <= xpd);
	FORI(n)
		FORL(ypd)
			y[ypd*i + l] = x[xpd*i + l];
}

#include "smapa.h"
SMART_PARAMETER(SRAND,0)
SMART_PARAMETER(PLAMBDA_CALC,0)

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

	float out[pdmax];
	int od = run_program_vectorially_at(out, p, x, 1, 1, pd, 0, 0);

	for (int i = 0; i < od; i++)
		fprintf(stderr, "result[%d] = %g\n", i, out[i]);

	return EXIT_SUCCESS;
}

int main_working(int c, char *v[])
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s in1 in2 ... \"plambda\"\n", *v);
		//                          0 1   2         c-1
		return EXIT_FAILURE;
	}

	struct plambda_program p[1];

	plambda_compile_program(p, v[c-1]);
	print_compiled_program(p);

	int n = c - 2;
	if (n != p->var->n)
		fail("the program expects %d variables but %d images "
					"were given", p->var->n, n);
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i+1], w + i, h + i, pd + i);
	FORI(n-1)
		if (w[0] != w[i+1] || h[0] != h[i+1])// || pd[0] != pd[i+1])
			fail("input images size mismatch");

	FORI(n)
		fprintf(stderr, "correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	srand(SRAND());

	int pdmax = PLAMBDA_MAX_PIXELDIM;
	float *out = xmalloc(*w * *h * pdmax * sizeof*out);
	//fprintf(stderr, "w h pd = %d %d %d\n", *w, *h, *pd);
	int opd = run_program_vectorially(out, pdmax, p, x, *w, *h, pd);
	//fprintf(stderr, "opd = %d\n", opd);

	float *rout = xmalloc(*w * *h *opd * sizeof*rout);
	shrink_components(rout, out, *w * *h, opd, pdmax);

	iio_save_image_float_vec("-", rout, *w, *h, opd);


	FORI(n) free(x[i]);
	free(out);
	free(rout);
	collection_of_varnames_end(p->var);

	return EXIT_SUCCESS;
}

int main_working2(int c, char *v[])
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s in1 in2 ... \"plambda\"\n", *v);
		//                          0 1   2         c-1
		return EXIT_FAILURE;
	}

	struct plambda_program p[1];

	plambda_compile_program(p, v[c-1]);
	print_compiled_program(p);

	int n = c - 2;
	if (n != p->var->n)
		fail("the program expects %d variables but %d images "
					"were given", p->var->n, n);
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i+1], w + i, h + i, pd + i);
	FORI(n-1)
		if (w[0] != w[i+1] || h[0] != h[i+1])// || pd[0] != pd[i+1])
			fail("input images size mismatch");

	FORI(n)
		fprintf(stderr, "correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	srand(SRAND());

	int pdreal = eval_dim(p, x, pd);

	//int pdmax = PLAMBDA_MAX_PIXELDIM;
	float *out = xmalloc(*w * *h * pdreal * sizeof*out);
	//fprintf(stderr, "w h pd = %d %d %d\n", *w, *h, *pd);
	int opd = run_program_vectorially(out, pdreal, p, x, *w, *h, pd);
	assert(opd == pdreal);
	//fprintf(stderr, "opd = %d\n", opd);

	iio_save_image_float_vec("-", out, *w, *h, opd);


	FORI(n) free(x[i]);
	free(out);
	collection_of_varnames_end(p->var);

	return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
	int (*f)(int c, char *v[]);
       	f = PLAMBDA_CALC()>0 ? main_calc : main_working2;
	return f(c,v);
}
