#include "fancy_image.h"
#include <stdio.h>
int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s img < points.txt\n", *v);
	struct fancy_image *a = fancy_image_open(v[1], "r");
	float x, y;
	while (2 == scanf("%g %g\n", &x, &y))
		for (int i = 0; i < a->pd; i++)
			printf("%g%c", fancy_image_getsample(a, x, y, i),
					i == a->pd - 1 ? '\n' : ' ');
	return 0;
}
