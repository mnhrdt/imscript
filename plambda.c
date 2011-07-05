#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"

#define PLAMBDA_MAX_TOKENS 0x100



#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif

#include "fragments.c"


#define PLAMBDA_CONSTANT 0
#define PLAMBDA_VARIABLE_SCALAR 1
#define PLAMBDA_VARIABLE_VECTOR 2
#define PLAMBDA_OPERATOR 3

struct predefined_function {
	void (*f)(void);
	char *name;
	int nargs;
	float value;
} global_table_of_predefined_functions[] = {
#define REGISTER_FUNCTION(x,n) {(void(*)(void))x, #x, n, 0}
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
#undef REGISTER_FUNCTION
	{NULL, "pi", 0, M_PI}
};


struct plambda_token {
	int type;
	int nargs;           // if type==operator, number of expected arguments
	float value;         // if type==constant, value
	int index;           // if type==variable, its index
	int component;       // if type==variable, index of selected component
	int displacement[2]; // if type==variable, relative displacement
};


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

// if the token is a valid identifier, return its length
// if the token is followed by modifiers, fill *endptr
static int token_is_identifier(const char *t, const char **endptr)
{
	*endptr = NULL;
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

struct predefined_function *identifier_is_predefined(const char *id)
{
	int n = sizeof(global_table_of_predefined_functions)/
		sizeof(global_table_of_predefined_functions[0]);
	struct predefined_function *r = global_table_of_predefined_functions;
	FORI(n)
		if (0 == strcmp(r[i].name, id))
			return r + i;
	return NULL;
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

static void identify_token(const char *tokke)
{
	char tok[1+strlen(tokke)];
	strncpy(tok, tokke, 1+strlen(tokke));

	if (false
			|| 0 == strcmp(tok, "+")
			|| 0 == strcmp(tok, "-")
			|| 0 == strcmp(tok, "/")
			|| 0 == strcmp(tok, "*")
			|| 0 == strcmp(tok, "^")
	   ) {
		fprintf(stderr, "TKN arithmetical operator '%c'\n", tok[0]);
		return;
	}

	float x;
	if (token_is_number(&x, tok))
	{
		fprintf(stderr, "TKN numeric constant: %g\n", x);
		return;
	}

	const char *endptr;
	if (token_is_identifier(tok, &endptr))
	{
		fprintf(stderr, "TKN identifier: \"%s\"\n", tok);
		struct predefined_function *p = identifier_is_predefined(tok);
		if (p)
			fprintf(stderr, "TKN predefined: %s (%d)\n",
					p->name, p->nargs);
		else {
			int comp, disp[2];
			parse_modifiers(endptr, &comp, disp, disp+1);
			fprintf(stderr, "TKN modifiers (%d,%d)[%d]\n",
					disp[0], disp[1], comp);
			fprintf(stderr, "TKN %s\n",
					comp<0 ? "VECTOR VARIABLE"
					: "SCALAR VARIABLE");
		}
		return;
	}

	fprintf(stderr, "TKN unrecognized\n");
}

static int plambda_tokenize(struct plambda_token *t, const char *str)
{
	// strdupa
	int slen = 1 + strlen(str);
	char s[slen];
	strncpy(s, str, slen);
	char *spacing = " ";

	int n = 0;
	char *tok = strtok(s, spacing);
	while (tok) {
		fprintf(stderr, "\nTOK[%d] = %s\n", n, tok);
		identify_token(tok);
		tok = strtok(NULL, spacing);
		n += 1;
	}
	return n;
}

int main(int c, char *v[])
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s in1 in2 ... \"plambda\"\n", *v);
		//                          0 1   2         c-1
		return EXIT_FAILURE;
	}
	int n = c - 2;
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i], w + i, h + i, pd + i);
	FORI(n-1)
		if (w[0] != w[i+1] || h[0] != h[i+1] || pd[0] != pd[i+1])
			error("input images size mismatch");

	struct plambda_token t[PLAMBDA_MAX_TOKENS];
	int ntok = plambda_tokenize(t, v[c-1]);

	FORI(n) free(x[i]);
	return EXIT_SUCCESS;
}
