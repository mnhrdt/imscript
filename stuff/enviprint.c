#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
	fprintf(stderr, "argc = %d\n", argc);
	for (int i = 0; i < argc; i++)
		fprintf(stderr, "argv[%d] = \"%s\"\n", i, argv[i]);
	for (char **s = environ; *s; s += 1)
		fprintf(stderr, "environ: \"%s\"\n", *s);
	return 33;
}
