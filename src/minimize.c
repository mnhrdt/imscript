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
// 	./minimize program "optvec" "startvec" "stepvec"
//
//

#define _POSIX_C_SOURCE 2

#include <stdio.h>
#include <stdlib.h>

#include "fail.c"


#include "minimize_gsl.c"


struct program {
	char *name;
	char *optset;
	int n;
};

static double run_program(struct program *p, double *x)
{
	static int run_idx = 0;

	int bufsize = 0x1000;
	char buf[bufsize];

	int n = snprintf(buf, bufsize, "%s %s", p->name, p->optset);
	for (int i = 0; i < p->n; i++)
		n += snprintf(buf+n, bufsize-n, " %f", x[i]);
	//snprintf(buf+n, bufsize-n, "");
	fprintf(stderr, "RUN[%d] %s", run_idx,  buf);

	FILE *f = popen(buf, "r");
	double y;
	if (1 != fscanf(f, "%lf\n", &y))
		fail("\ncommand [%s] did not print a number\n", buf);
	pclose(f);

	fprintf(stderr, "\t=> %g\n", y);

	run_idx += 1;
	return y;
}

static double run_program_as_of(double *x, int n, void *data)
{
	struct program *p = data;
	assert(p->n ==  n);
	return run_program(p, x);
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
	if (c != 5) {
		fprintf(stderr, "usage:\n\t"
		"%s program \"optvec\" \"parvec\" \"stepvec\"\n", *v);
	//       0  1         2          3          4
		return EXIT_FAILURE;
	}
	struct program p[1] = {{.name = v[1], .optset = v[2]}};
	int nmax = 100;
	double startvec[nmax], stepvec[nmax];
	int n = parse_doubles(startvec, nmax, v[3]);
	int ns = parse_doubles(stepvec, nmax, v[4]);
	if (n != ns)
		fail("dimensions %d and %d do not match", n, ns);
	p->n = n;

	//double x = run_program(p, startvec);
	//fprintf(stderr, "got %g\n", x);

	double result[n];
	int r = minimize_objective_function(result, startvec, stepvec,
			run_program_as_of, n, p);

	fprintf(stderr, "minimization returned %d\n", r);

	return EXIT_SUCCESS;
}
