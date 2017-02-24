#include <string.h>
int main(int c, char *v[])
{
#define MAIN(x) int main_ ## x(int, char**);\
	if (c>1 && !strcmp(v[1], #x)) return main_ ## x(c-1, v+1)
#include "all_mains.inc"
	return 1;
}
