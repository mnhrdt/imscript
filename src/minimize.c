// A program for minimizing other programs
//
//
//
// 1. CONTEXT
//
// Suppose that you have a program (or a script file) that is called
// in the following way:
//
// 	./program opt1 opt2 ... optM par1 par2 ... parN
//
// where par1 ... parN are floating-point numbers and opt1 ... optN are
// arbitrary strings.
//
// Suppose that running this program prints a single floating point
// number to stdout, which is computed deterministically from its
// arguments.  This number is called the "result" of the program.
//
// Suppose that you want to find the set of values of par1 ... parN
// where the result attains its minimum (for a fixed choice of opt1
// ... optM).
//
// Then, this program is for you.  Although global minimization is
// unsolvable in general, doing it automatically is more confortable
// than doing it by hand, and faster than by exhaustive search on a
// fixed grid.
//
//
// 2. SYNTAX
//
// 	./minimize program "optset" "parset0" "parset1" ... "parsetN"
//
//

#define _POSIX_C_SOURCE 2

#include <stdio.h>
#include <stdlib.h>

#include "fail.c"

struct program {
	char *name;
	char *optset;
	int n;
};

static double run_program(struct program *p, double *x)
{
	int bufsize = 0x1000;
	char buf[bufsize];

	int n = snprintf(buf, bufsize, "%s \"%s\" \"", p->name, p->optset);
	for (int i = 0; i < p->n; i++)
		n += snprintf(buf+n, bufsize-n, " %f", x[i]);
	snprintf(buf+n, bufsize-n, "\"");
	fprintf(stderr, "RUN %s\n", buf);

	FILE *f = popen(buf, "r");
	double y;
	if (1 != fscanf(f, "%lf\n", &y))
		fail("command [%s] did not print a number\n", buf);
	pclose(f);

	return y;
}

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

int main(int c, char *v[])
{
	if (c < 4) {
		fprintf(stderr, "usage:\n\t"
		"%s program \"optset\" \"parset0\" ... \"parsetN\"\n", *v);
	//       0  1         2          3               c-1
		return EXIT_FAILURE;
	}
	int nmax = c  - 3;
	struct program p = {.name = v[1], .optset = v[2], .n = nmax};
	for (int i = 3; i < c; i++) {
		double t[nmax];
		int n = parse_doubles(t, nmax, v[i]);
		if (n != nmax)
			fail("you must supply %d vectors of length %d,\n"
			"however, your %dth vector has length %d\n",
						nmax+1, nmax, i-2, n);
		double x = run_program(&p, t);
		fprintf(stderr, "returned x = %g\n", x);
	}

	return EXIT_SUCCESS;
}
