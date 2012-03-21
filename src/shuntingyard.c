#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fail.c"
#include "xmalloc.c"


enum token_id {
	TOK_VALUE, // constant, scalar, vector
	TOK_FUNCTION,
	TOK_COMMA,
	TOK_OPERATOR,
	TOK_LEFTPAR,
	TOK_RIGHTPAR,
	TOK_SEMICOLON,
	TOK_ASSIGN,
};

static char *token_id_name(enum token_id t)
{
	switch(t) {
#define X(s) case s: return #s
	X(TOK_VALUE);
	X(TOK_FUNCTION);
	X(TOK_COMMA);
	X(TOK_OPERATOR);
	X(TOK_SEMICOLON);
	X(TOK_LEFTPAR);
	X(TOK_RIGHTPAR);
	X(TOK_ASSIGN);
#undef X
	default: fail("unrecognized token id %d", t);
	}
}



#define MAX_TOK_LEN 33
struct token {
	enum token_id type;

	char s[MAX_TOK_LEN];

};


static enum token_id get_token_type(char *str)
{
	if (0 == strcmp(str, "(")) return TOK_LEFTPAR;
	if (0 == strcmp(str, ")")) return TOK_RIGHTPAR;
	if (0 == strcmp(str, ",")) return TOK_COMMA;
	if (0 == strcmp(str, ";")) return TOK_SEMICOLON;
	if (0 == strcmp(str, ":=")) return TOK_ASSIGN;
	if (0 == strcmp(str, "+")) return TOK_OPERATOR;
	if (0 == strcmp(str, "-")) return TOK_OPERATOR;
	if (0 == strcmp(str, "*")) return TOK_OPERATOR;
	if (0 == strcmp(str, "/")) return TOK_OPERATOR;
	if (0 == strcmp(str, "^")) return TOK_OPERATOR;
	if (0 == strcmp(str, "sin")) return TOK_FUNCTION;
	if (0 == strcmp(str, "cos")) return TOK_FUNCTION;
	if (0 == strcmp(str, "exp")) return TOK_FUNCTION;
	if (0 == strcmp(str, "log")) return TOK_FUNCTION;
	if (0 == strcmp(str, "atan2")) return TOK_FUNCTION;
	if (0 == strcmp(str, "pow")) return TOK_FUNCTION;
	if (0 == strcmp(str, "hypot")) return TOK_FUNCTION;
	return TOK_VALUE;
}

static int precedence(struct token *t)
{
	assert(t->type == TOK_OPERATOR);
	switch(t->s[0]) {
	case '+': return 1;
	case '-': return 1;
	case '*': return 3;
	case '/': return 3;
	case '^': return 5;
	default: assert(false);
	}
}

static int tokenize(struct token *t, char *str, int nmax)
{
	char s[1+strlen(str)];
	strcpy(s, str);
	char *spacing = " \n\t";

	int n = 0;
	char *tok = strtok(s, spacing);
	while(tok) {
		t[n].type = get_token_type(tok);
		strncpy(t[n].s, tok, MAX_TOK_LEN);
		fprintf(stderr, "TOK \"%s\" has type %s\n",
				t[n].s,
				token_id_name(t[n].type));
		tok = strtok(NULL, spacing);
		n += 1;
	}
	return n;
}

struct stack_of_operators {
	int n;
	struct token t[0];
};

static void pstack(struct stack_of_operators *s)
{
	fprintf(stderr, "stack:");
	for (int i = 0; i < s->n; i++)
		fprintf(stderr, " \"%s\"", s->t[i].s);
	fprintf(stderr, "\n");
}

static void push(struct stack_of_operators *s, struct token *t)
{
	s->t[s->n] = *t;
	s->n += 1;
}

static struct token *top(struct stack_of_operators *s)
{
	if (s->n <= 0)
		return NULL;
	else
		return s->t + s->n - 1;
}

static struct token *pop(struct stack_of_operators *s)
{
	struct token *r = top(s);
	if (r)
		s->n -= 1;
	return r;
}

static void queue(struct token *t)
{
	fprintf(stderr, "QUEUE \"%s\"\n", t->s);
}

static void shunting_yard(char *s)
{
	int ns = strlen(s);
	struct token *t = xmalloc(ns * sizeof*t), *tt;

	int nt = tokenize(t, s, ns);

	struct stack_of_operators *stack = xmalloc(sizeof*stack + ns*sizeof*t);
	stack->n = 0;

	for (int i = 0; i < nt; i++) {
		struct token *tok = t + i;
		//fprintf(stderr, "processing token \"%s\"\n", tok->s);
		//pstack(stack);
		switch(tok->type) {
		case TOK_VALUE:
			queue(tok);
			break;
		case TOK_FUNCTION:
			push(stack, tok);
			break;
		case TOK_COMMA:
			while ((tt = top(stack)) && tt->type != TOK_LEFTPAR)
				queue(pop(stack));
			if (!tt) fail("mismatched parentheses\n");
			break;
		case TOK_OPERATOR:
			while ((tt = top(stack)) && tt->type == TOK_OPERATOR
					&& precedence(tok) <= precedence(tt))
				queue(pop(stack));
			push(stack, tok);
			break;
		case TOK_LEFTPAR:
			push(stack, tok);
			break;
		case TOK_RIGHTPAR:
			while ((tt = top(stack)) && tt->type != TOK_LEFTPAR)
				queue(pop(stack));
			if (!tt) fail("mismatched parens\n");
			assert(tt->type == TOK_LEFTPAR);
			pop(stack);
			if ((tt = top(stack)) && tt->type == TOK_FUNCTION)
				queue(pop(stack));
			break;
		default:
			fail("don't know how to deal with %s\n",
					token_id_name(tok->type));
		}
	}

	//fprintf(stderr, "no more tokens\n");
	//pstack(stack);
	while (tt = pop(stack))
		if (tt->type == TOK_LEFTPAR || tt->type == TOK_RIGHTPAR)
			fail("unmatched parentheses\n");
		else
			queue(tt);

	xfree(t);
}

int main(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s expression\n", *v);
		return EXIT_FAILURE;
	}

	shunting_yard(v[1]);

	return EXIT_SUCCESS;
}
