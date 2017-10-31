#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t"
			"%s threshold < in.txt > out.txt\n", *v);
	//                0 1
	float t = atof(v[1]);

	float x[4];
	while (4 == scanf("%g %g %g %g\n", x+0, x+1, x+2, x+3))
		if (hypot(x[0] - x[2], x[1] - x[3]) > t)
			printf("%g %g %g %g\n", x[0], x[1], x[2], x[3]);

	return 0;
}
