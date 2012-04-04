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
	if (0 == strcmp(str, ":")) return TOK_ASSIGN;
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


static void add_to_token(char *tok, int c)
{
	int n = strlen(tok);
	tok[n] = c;
}

static void emit_token(char *tok)
{
	fprintf(stderr, "TOKEN \"tok\"\n", tok);
	for (int i = 0; i < MAX_TOK_LEN; i++)
		tok[i] = 0;
}

// hand-written naive lexer
static void experimental_tokenization(char *input_string)
{
	// character types
	enum {
		T_UNKNOWN,
		T_ALNUM,      // letters and numbers
		T_WHITESPACE, // space, tab, newline
		T_PUNCT,      // single-character punctuation and operators
		T_MINUS,      // -
		T_PLUS,       // +
		T_OPENCLY,    // {
		T_CLOSECLY,   // }
		T_OPENBRA,    // [
		T_CLOSEBRA,   // ]
		T_COMMA,      // ,
		T_DOT,        // .
		T_COLON,      // :
		T_END,	// '\0'
		T_N
	};


	// fill table of character types
	int t[0x100];
	for (int i = 0; i < 0x100; i++)  t[i] = T_UNKNOWN;
	for (int i = 'a'; i <= 'z'; i++) t[i] = t[toupper(i)] = T_ALNUM;
	for (int i = '0'; i <= '9'; i++) t[i] = T_ALNUM;
	t['+'] = T_PLUS;
	t['-'] = T_MINUS;
	t['{'] = T_OPENCLY;
	t['}'] = T_CLOSECLY;
	t['['] = T_OPENBRA;
	t[']'] = T_CLOSEBRA;
	t[','] = T_COMMA;
	t['.'] = T_DOT;
	t[':'] = T_COLON;
	t['<']=t['>']=t['*']=t['(']=t[')']=t['/']=t['^'] = T_PUNCT;
	t['\0'] = T_END;


	// state machine states (codifying the state at the last read symbol)
	enum {
		S_START,
		S_IDENT, // inside constant, identifier or its modifiers
		S_PUNCT, // punctuation and single-letter operators
		S_SPACE, // whitespace
		S_END,
		S_N
       	};

	// state machine actions
	enum {
		A_NOP,
		A_EMIT,
		A_FAIL
	}

	// build state machine: default actions and transitions
	int fsm[T_N][S_N][2];
	for (int j = 0; j < T_N; j++)
	for (int i = 0; i < S_N; i++) {
		fsm[j][i][0] = A_FAIL;
		fsm[j][i][1] = S_END;
	}

	// build state machine


	// run state machine
	int state = LEX_IN_START;
	char tok[MAX_TOK_LEN] = {0};

	int c;
	while ((c = *input_string++)) {
		assert(c >= 0 && c < 0x100);
		assert(state != S_END);
		fprintf(stderr, "read character '%c'\n", c);

		int tc = t[c];
		int action = fsm[tc][state][0];
		switch(action) {
		case A_NOP: break;
		case A_EMIT: emit_token(tok); break;
		case A_FAIL: fail("caca\n"); break;
		default: fail("merda\n");
		}
		state = fsm[tc][state][0];
		add_to_token(tok, c);
	}

}

static int tokenize(struct token *t, char *str, int nmax)
{
	experimental_tokenization(str);
	exit(2);
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
