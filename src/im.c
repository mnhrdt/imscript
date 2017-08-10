#include <string.h>
int main(int c, char *v[])
{
#define MAIN(x) int main_ ## x(int, char**);\
	if (c>1 && !strcmp(v[1], #x)) return main_ ## x(c-1, v+1)
#include "all_mains.inc"
//	MAIN(plambda);
//	MAIN(qeasy);
//	MAIN(blur);
//	MAIN(vecov);
//	MAIN(nnint);
//	MAIN(srmatch);
//	MAIN(ransac);
//	MAIN(iion);
	return 1;
}
