#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "xfopen.c"

static int getlinen(char *l, int n, FILE *f)
{
	int c, i = 0;
	while (i < n-1 && (c = fgetc(f)) != EOF && c != '\n')
		l[i++] = c;
	l[i] = '\0';
	return i;
}

static int insideP(int w, int h, int x, int y)
{
	return x >= 0 && y >= 0 && x < w && y < h;
}

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
				"%s hA.tiff hB.tiff [pairs2d [pairs3d]]\n", *v);
		//                0 1       2        3        4
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *filename_in  = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int wa, ha, wb, hb;
	float *a = iio_read_image_float(filename_a, &wa, &ha);
	float *b = iio_read_image_float(filename_b, &wb, &hb);
	FILE *fi = xfopen(filename_in, "r");
	FILE *fo = xfopen(filename_out, "w");

	int n, lmax = 10000;
	char line[lmax];
	while (n = getlinen(line, lmax, fi))
	{
		double xy[4], xyz[6];
		int r = sscanf(line, "%lf %lf %lf %lf", xy, xy+1, xy+2, xy+3);
		if (r != 4) continue;
		int ia = lrint(xy[0]);
		int ja = lrint(xy[1]);
		int ib = lrint(xy[2]);
		int jb = lrint(xy[3]);
		if (!insideP(wa, ha, ia, ja)) continue;
		if (!insideP(wb, hb, ib, jb)) continue;
		double va = a[wa * ja + ia];
		double vb = b[wb * jb + ib];
		if (!isfinite(va)) continue;
		if (!isfinite(vb)) continue;
		fprintf(fo, "%lf %lf %lf %lf %lf %lf\n",
				xy[0], xy[1], va, xy[2], xy[3], vb);
	}

	free(a);
	free(b);
	xfclose(fo);
	xfclose(fi);
	return 0;
}
