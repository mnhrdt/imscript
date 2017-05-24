#include <stdio.h>
#include <string.h>
#define MAIN(x) int main_ ## x(int, char**);
#include "all_mains.inc"
#undef MAIN
static struct {
	char *n;
	int (*f)(int,char **);
} t[] = {
#define MAIN(x) { #x, main_ ## x },
#include "all_mains.inc"
#undef MAIN
	{0, 0}
};
int main(int c, char *v[])
{
	int i = 0;
	while (t[i].n)
	{
		if (c == 1)
		{
			char *newv[] = {t[i].n, "--help-oneliner", 0};
			printf("%s\n", t[i].n);
			//fflush(stdout);
			//t[i].f(2, newv);
		}
		if (c > 1 && !strcmp(v[1], t[i].n))
			return t[i].f(c - 1, v + 1);
		i = i + 1;
	}
	return i;
}
