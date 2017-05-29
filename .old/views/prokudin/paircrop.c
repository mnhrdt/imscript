#include <stdio.h>
#include <stdlib.h>

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s x0 y0 xf yf < pairs\n", *v);
		//                          0 1  2  3  4
		return 1;
	}
	int x0 = atoi(v[1]);
	int y0 = atoi(v[2]);
	int xf = atoi(v[3]);
	int yf = atoi(v[4]);

	float k[4];
	while (4 == scanf("%g %g %g %g\n", k, k+1, k+2, k+3))
		if (x0 <= k[0] && k[0] <= xf    &&    y0 <= k[1] && k[1] <= yf)
			printf("%g %g %g %g\n", k[0], k[1], k[2], k[3]);

	return 0;
}
