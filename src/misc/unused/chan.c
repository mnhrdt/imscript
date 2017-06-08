#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#define xmalloc malloc

// chan "list of desired channels" [in [out]]
//
// Available channels:
// R red
// G green
// B blue
// r R/(R+G+B)
// g G/(R+G+B)
// b B/(R+G+B)
// H hue
// S saturation
// L lightness
// l luminance
// V value
// I intensity
// A (from LAB)
// B (from LAB)
// X
// Y
// pca1
// pca2
// pca3
// tex1
// tex2
// tex3
// ...
int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s idx [in [out]]\n", *v);
		//                         0  1    2   3
		return EXIT_FAILURE;
	}
	char *in = c > 2 ? v[2] : "-";
	char *out = c > 3 ? v[3] : "-";
	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	struct component_list l[1];
	int nc = parse_components(l, v[1]);
	float *y = xmalloc(w * h * nc * sizeof*y);
	fill_requested_components(y, l, x);
	iio_write_image_float(out, y, w, h);
	return EXIT_SUCCESS;
}
