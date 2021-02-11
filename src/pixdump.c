#include <stdio.h> // fprintf, scanf, printf
#include "iio.h"
int main(int c, char *v[])
{
	if (c > 2 || (c==2 && v[1][0]=='-'))
		return fprintf(stderr, "usage:\n echo x y | %s img\n", *v);

	int w, h, d;
	float *x = iio_read_image_float_vec(v[1], &w, &h, &d);
	int p, q;
	while (2 == scanf("%d %d\n", &p, &q))
		if (p >= 0 && p < w && q >= 0 && q < h)
			for (int i = 0; i < d; i++)
				printf("%g%c", x[(q*w+p)*d+i], i<d-1?' ':'\n');
	return 0;
}
