// NAME
// 	plambda - a RPN calculator for image pixels
//
// SYNOPSIS
// 	plambda a.pnb b.png c.png ... "lambda expression" > output.png
//
// DESCRIPTION
// 	Plambda applies an expression to all the pixels of a collection of
// 	images, and produces a single output image.  Each input image
// 	corresponts to one of the variables of the expression (in alphabetical
// 	order).  There are modifiers to the variables that allow access to the
// 	values of neighboring pixels.
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
//
//
// TODO: complete support for color images (incomplete, yet not called)
// TODO: implement support for nd images (easy, mostly fix i/o problems)
// TODO: implement shunting-yard algorithm to admit infix notation


#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#define PLAMBDA_MAX_TOKENS 0x100
#define PLAMBDA_MAX_VARLEN 0x100
#define PLAMBDA_MAX_PIXELDIM 40



#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif
#ifndef FORJ
#define FORJ(n) for(int j=0;j<(n);j++)
#endif
#ifndef FORK
#define FORK(n) for(int k=0;k<(n);k++)
#endif

#include "fragments.c"


#define PLAMBDA_CONSTANT 0   // numeric constant
#define PLAMBDA_SCALAR 1     // pixel component
#define PLAMBDA_VECTOR 2     // whole pixel
#define PLAMBDA_OPERATOR 3   // function
#define PLAMBDA_COLONVAR 4   // colon-type variable

static double sum_two_doubles      (double a, double b) { return a + b; }
static double substract_two_doubles(double a, double b) { return a - b; }
static double multiply_two_doubles (double a, double b) { return a * b; }
static double divide_two_doubles   (double a, double b) { return a / b; }

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
	REGISTER_FUNCTIONN(pow,"^",2),
	REGISTER_FUNCTIONN(sum_two_doubles,"+",2),
	REGISTER_FUNCTIONN(divide_two_doubles,"/",2),
	REGISTER_FUNCTIONN(multiply_two_doubles,"*",2),
	REGISTER_FUNCTIONN(substract_two_doubles,"-",2),
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

	char *tmphack;
};

struct collection_of_varnames {
	int n;
	char *t[PLAMBDA_MAX_TOKENS];
	//int count[PLAMBDA_MAX_TOKENS];
};

struct plambda_program {
	int n;
	struct plambda_token t[PLAMBDA_MAX_TOKENS];
	struct collection_of_varnames var[1];
};


static float apply_function(struct predefined_function *f, float *v)
{
	switch(f->nargs) {
	case 0: return f->value;
	case 1: return ((double(*)(double))(f->f))(v[0]);
	case 2: return ((double(*)(double,double))f->f)(v[1], v[0]);
	case 3: error("bizarre");
	}
	return 0;
}

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
	default: error("unrecognized colonvar \":%c\"", c);
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

static int token_is_colonvar(const char *t)
{
	if (t[0] != ':') return 0;
	if (isalpha(t[1]) && t[2]=='\0') return t[1];
	return 0;
}

// if the token is a valid word, return its length
// if the token is followed by modifiers, fill *endptr
static int token_is_word(const char *t, const char **endptr)
{
	*endptr = NULL;
	if ((*t=='+'||*t=='-'||*t=='/'||*t=='^'||*t=='*')&&t[1]=='\0')
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
			error("caca");
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
static void process_token(struct plambda_program *p, const char *tokke)
{
	char tok[1+strlen(tokke)];             // the string of the token
	strcpy(tok, tokke);
	struct plambda_token *t = p->t + p->n; // the compiled token

	float x;
	if (token_is_number(&x, tok))
	{
		//fprintf(stderr, "TKN numeric constant: %g\n", x);
		t->type = PLAMBDA_CONSTANT;
		t->value = x;
		goto endtok;
	}

	int colonvar = token_is_colonvar(tok);
	if (colonvar)
	{
		t->type = PLAMBDA_COLONVAR;
		t->colonvar = colonvar;
		goto endtok;
	}

	const char *endptr;
	if (token_is_word(tok, &endptr))
	{
		//fprintf(stderr, "TKN word: \"%s\"\n", tok);
		int idx = word_is_predefined(tok);
		if (idx < 0) {
			char varname[PLAMBDA_MAX_VARLEN+1];
			int varlen = endptr-tok;
			if (varlen >= PLAMBDA_MAX_VARLEN)
				varlen = PLAMBDA_MAX_VARLEN;
			FORI(varlen) varname[i] = tok[i];
			varname[varlen] = '\0';
			int comp, disp[2];
			//fprintf(stderr, "TKN varname = %s\n", varname);
			t->tmphack =collection_of_varnames_add(p->var, varname);
			parse_modifiers(endptr, &comp, disp, disp+1);
			//fprintf(stderr, "TKN modifiers (%d,%d)[%d]\n",
			//		disp[0], disp[1], comp);
			//fprintf(stderr, "TKN %s\n",
			//		comp<0 ? "VECTOR VARIABLE"
			//		: "SCALAR VARIABLE");
			t->type = comp<0 ? PLAMBDA_VECTOR : PLAMBDA_SCALAR;
			t->component = comp;
			t->displacement[0] = disp[0];
			t->displacement[1] = disp[1];
		} else {
			struct predefined_function *f =
				global_table_of_predefined_functions + idx;
			//fprintf(stderr, "TKN predefined function: %s (%d)\n",
			//		f->name, f->nargs);
			t->type = PLAMBDA_OPERATOR;
			t->index = idx;
		}
		goto endtok;
	}

endtok:
	p->n += 1;
}

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
				error("unexpected bad variable \"%s\"",
								t->tmphack);
		}
	}
}

//static void identify_token(struct collection_of_varnames *cv, const char *tokke)
//{
//	char tok[1+strlen(tokke)];
//	strncpy(tok, tokke, 1+strlen(tokke));
//
//	if (false
//			|| 0 == strcmp(tok, "+")
//			|| 0 == strcmp(tok, "-")
//			|| 0 == strcmp(tok, "/")
//			|| 0 == strcmp(tok, "*")
//			|| 0 == strcmp(tok, "^")
//	   ) {
//		fprintf(stderr, "TKN arithmetical operator '%c'\n", tok[0]);
//		return;
//	}
//
//	float x;
//	if (token_is_number(&x, tok))
//	{
//		fprintf(stderr, "TKN numeric constant: %g\n", x);
//		return;
//	}
//
//	const char *endptr;
//	if (token_is_word(tok, &endptr))
//	{
//		fprintf(stderr, "TKN word: \"%s\"\n", tok);
//		struct predefined_function *p = word_is_predefined(tok);
//		if (p)
//			fprintf(stderr, "TKN predefined function: %s (%d)\n",
//					p->name, p->nargs);
//		else {
//			char varname[PLAMBDA_MAX_VARLEN+1];
//			int varlen = endptr-tok;
//			if (varlen >= PLAMBDA_MAX_VARLEN)
//				varlen = PLAMBDA_MAX_VARLEN;
//			FORI(varlen) varname[i] = tok[i];
//			varname[varlen] = '\0';
//			int comp, disp[2];
//			fprintf(stderr, "TKN varname = %s\n", varname);
//			collection_of_varnames_add(cv, varname);
//			parse_modifiers(endptr, &comp, disp, disp+1);
//			fprintf(stderr, "TKN modifiers (%d,%d)[%d]\n",
//					disp[0], disp[1], comp);
//			fprintf(stderr, "TKN %s\n",
//					comp<0 ? "VECTOR VARIABLE"
//					: "SCALAR VARIABLE");
//		}
//		return;
//	}
//
//	fprintf(stderr, "TKN unrecognized\n");
//}

//static int plambda_tokenize(struct plambda_token *t, const char *str)
//{
//	// strdupa
//	int slen = 1 + strlen(str);
//	char s[slen];
//	strncpy(s, str, slen);
//	char *spacing = " ";
//
//	struct collection_of_varnames cvars[1];
//	collection_of_varnames_init(cvars);
//
//	int n = 0;
//	char *tok = strtok(s, spacing);
//	while (tok) {
//		fprintf(stderr, "\nTOK[%d] = %s\n", n, tok);
//		identify_token(cvars, tok);
//		tok = strtok(NULL, spacing);
//		n += 1;
//	}
//
//	collection_of_varnames_sort(cvars);
//	fprintf(stderr, "\nwe got %d variables:\n", cvars->n);
//	FORI(cvars->n)
//		fprintf(stderr, "\tVAR[%d] = \"%s\"\n", i, cvars->t[i]);
//
//	collection_of_varnames_end(cvars);
//
//	return n;
//}

static void plambda_compile_program(struct plambda_program *p, const char *str)
{
	char s[1+strlen(str)];
	strcpy(s, str);
	char *spacing = " ";

	collection_of_varnames_init(p->var);
	p->n = 0;
	int n = 0;
	char *tok = strtok(s, spacing);
	while (tok) {
		//fprintf(stderr, "\nTOK[%d] = %s\n", n, tok);
		process_token(p, tok);
		tok = strtok(NULL, spacing);
		n += 1;
	}

	collection_of_varnames_sort(p->var);
	//fprintf(stderr, "\nwe got %d variables:\n", p->var->n);
	//FORI(p->var->n)
	//	fprintf(stderr, "\tVAR[%d] = \"%s\"\n", i, p->var->t[i]);

	unhack_varnames(p);
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
			fprintf(stderr, "operator %s", f->name);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "The program uses %d variables:\n", p->var->n);
	FORI(p->var->n)
		fprintf(stderr, "VARIABLE[%d] = \"%s\"\n", i, p->var->t[i]);
}

struct value_stack {
	int n;
	float t[PLAMBDA_MAX_TOKENS];
};

struct value_vstack {
	int n;
	int d[PLAMBDA_MAX_TOKENS];
	float t[PLAMBDA_MAX_TOKENS][PLAMBDA_MAX_PIXELDIM];
};

static float stack_pop(struct value_stack *s)
{
	if (s->n > 0) {
		s->n -= 1;
		return s->t[s->n];
	} else error("popping from empty stack");
}

static int vstack_pop_vector(float *val, struct value_vstack *s)
{
	if (s->n > 0) {
		s->n -= 1;
		int d = s->d[s->n];
		FORI(d) val[i] = s->t[s->n][i];
		return d;
	} else error("popping from empty stack");
}

static float vstack_pop_scalar(struct value_vstack *s)
{
	if (s->n > 0) {
		s->n -= 1;
		if (s->d[s->n] != 1)
			error("trying to pop a scalar but found a vector");
		return s->t[s->n][0];
	} else error("popping from empty stack");
}

static void stack_push(struct value_stack *s, float x)
{
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		s->t[s->n] = x;
		s->n += 1;
	} else error("full stack");
}

static void vstack_push_vector(struct value_vstack *s, float *v, int n)
{
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		FORI(n)
			s->t[s->n][i] = v[i];
		s->n += 1;
	} else error("full stack");
}

static void vstack_push_scalar(struct value_vstack *s, float *v, int n)
{
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		FORI(n)
			s->t[s->n][i] = v[i];
		s->n += 1;
	} else error("full stack");
}

static void stack_print(FILE *f, struct value_stack *s)
{
	FORI(s->n)
		fprintf(f, "STACK[%d/%d]: %g\n", 1+i, s->n, s->t[i]);
}

static void vstack_print(FILE *f, struct value_vstack *s)
{
	FORI(s->n) {
		fprintf(f, "STACK[%d/%d]: {%d}", 1+i, s->n, s->d[i]);
		FORJ(s->d[i])
			fprintf(f, " %g", s->t[i][j]);
		fprintf(f, "\n");
	}
}


#include "getpixel.c"

static float run_program_at_pixel(struct plambda_program *p,
		float **val, int w, int h, int ai, int aj)
{
	struct value_stack s[1];
	s->n = 0;
	FORI(p->n) {
		struct plambda_token *t = p->t + i;
		switch(t->type) {
		case PLAMBDA_CONSTANT:
			stack_push(s, t->value);
			break;
		case PLAMBDA_COLONVAR: {
			float x = eval_colonvar(w, h, ai, aj, t->colonvar);
			stack_push(s, x);
			break;
				       }
		case PLAMBDA_SCALAR:
		case PLAMBDA_VECTOR: {
			float *img = val[t->index];
			int dai = ai + t->displacement[0];
			int daj = aj + t->displacement[1];
			float x = getsample_0(img, w, h, 1, dai, daj, 0);
			stack_push(s, x);
			//stack_push(s, val[t->index]);
				     }
			break;
		case PLAMBDA_OPERATOR: {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			float param[f->nargs];
			FORJ(f->nargs)
				param[j] = stack_pop(s);
			float x = apply_function(f, param);
			stack_push(s, x);
				       }
			break;
		default:
			error("unknown tag type %d", t->type);
		}
	}
	return stack_pop(s);
}

static void run_program(float *out,
		struct plambda_program *p, float **val, int w, int h)
{
	FORJ(h) FORI(w) {
//static float run_program_at_pixel(struct plambda_program *p,
//		float **val, int w, int h, int ai, int aj)
		float r = run_program_at_pixel(p, val, w, h, i, j);
		setsample_0(out, w, h, 1, i, j, 0, r);
	}
}

static void run_program_numbers(struct plambda_program *p, float *val)
{
	struct value_stack s[1];
	s->n = 0;
	FORI(p->n) {
		fprintf(stderr, "\n\nstack before token %d (%d):\n", i, s->n);
		stack_print(stderr, s);
		struct plambda_token *t = p->t + i;
		switch(t->type) {
		case PLAMBDA_CONSTANT:
			stack_push(s, t->value);
			break;
		case PLAMBDA_SCALAR:
		case PLAMBDA_VECTOR:
			stack_push(s, val[t->index]);
			break;
		case PLAMBDA_OPERATOR: {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			float param[f->nargs];
			FORJ(f->nargs)
				param[j] = stack_pop(s);
			float x = apply_function(f, param);
			stack_push(s, x);
				       }
			break;
		default:
			error("unknown tag type %d", t->type);
		}
	}
	fprintf(stderr, "\n\nstack at the end (%d elements):\n", s->n);
	stack_print(stderr, s);
}

int main(int c, char *v[])
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
		error("the program expects %d variables but %d images "
					"were given", p->var->n, n);
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i+1], w + i, h + i, pd + i);
	FORI(n-1)
		if (pd[i] != 1)
			error("only gray images supported so far");
	FORI(n-1)
		if (w[0] != w[i+1] || h[0] != h[i+1] || pd[0] != pd[i+1])
			error("input images size mismatch");

	FORI(n)
		fprintf(stderr, "correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	float *out = xmalloc(*w * *h * 1 * sizeof*out);
	run_program(out, p, x, *w, *h);
	iio_save_image_float_vec("-", out, *w, *h, 1);


	//FORI(n) free(x[i]);
	collection_of_varnames_end(p->var);

	return EXIT_SUCCESS;
}
