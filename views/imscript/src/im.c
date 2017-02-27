#include <stdio.h>
#include <string.h>
#define MAIN(x) int main_ ## x(int, char**);
#include "all_mains.inc"
static struct {
	char *n;
	int (*f)(int,char **);
} t[] = {
#undef MAIN
#define MAIN(x) { #x, main_ ## x },
#include "all_mains.inc"
	{0, 0}
};
int main(int c, char *v[])
{
	int i = 0;
	while (t[i].n)
	{
		if (c == 1)
			fprintf(stderr, "%s\n", t[i].n);
		if (c > 1 && !strcmp(v[1], t[i].n))
			return t[i].f(c - 1, v + 1);
		i = i + 1;
	}
	return i;
}
