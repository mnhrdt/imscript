#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#include "xfopen.c"
#include "parsenumbers.c"

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

static void set_identity(double H[9])
{
	H[0] = H[4] = H[8] = 1;
	H[1] = H[2] = H[5] = 0;
	H[3] = H[6] = H[7] = 0;
}

// y = H(x)
static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

// compute the inverse homography (inverse of a 3x3 matrix)
static double invert_homography(double invH[3][3], double H[3][3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = ( a[4] * a[8] - a[5] * a[7] ) / det;
	r[1] = ( a[2] * a[7] - a[1] * a[8] ) / det;
	r[2] = ( a[1] * a[5] - a[2] * a[4] ) / det;
	r[3] = ( a[5] * a[6] - a[3] * a[8] ) / det;
	r[4] = ( a[0] * a[8] - a[2] * a[6] ) / det;
	r[5] = ( a[2] * a[3] - a[0] * a[5] ) / det;
	r[6] = ( a[3] * a[7] - a[4] * a[6] ) / det;
	r[7] = ( a[1] * a[6] - a[0] * a[7] ) / det;
	r[8] = ( a[0] * a[4] - a[1] * a[3] ) / det;
	return det;
}

int main(int c, char *v[])
{
	if (c != 5 && c != 6 && c != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s Axyz.tiff Bxyx.tiff homA homB [pairs2d [pairs3d]]\n",*v);
		//0 1         2         3    4     5        6
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *homstring_a = v[3];
	char *homstring_b = v[4];
	char *filename_in  = c > 5 ? v[5] : "-";
	char *filename_out = c > 6 ? v[6] : "-";

	double Ha[9], Hb[9], iHa[3][3], iHb[3][3];
	int nHa = read_n_doubles_from_string(Ha, homstring_a, 9);
	int nHb = read_n_doubles_from_string(Hb, homstring_b, 9);
	if (nHa != 9) set_identity(Ha);
	if (nHb != 9) set_identity(Hb);
	invert_homography(iHa, (void*)Ha);
	invert_homography(iHb, (void*)Hb);

	int wa, ha, pda, wb, hb, pdb;
	float *a = iio_read_image_float_vec(filename_a, &wa, &ha, &pda);
	float *b = iio_read_image_float_vec(filename_b, &wb, &hb, &pdb);
	FILE *fi = xfopen(filename_in, "r");
	FILE *fo = xfopen(filename_out, "w");
	if (pda != pdb)
		return fprintf(stderr, "input images dimension mismatch\n");
	if (pda != 1 && pda != 3)
		return fprintf(stderr, "input images should be h or xyz\n");

	int n, lmax = 10000;
	char line[lmax];
	while ((n = getlinen(line, lmax, fi)))
	{
		double m[4], p[2], q[2];
		int r = sscanf(line, "%lf %lf %lf %lf", m, m + 1, m + 2, m + 3);
		if (r != 4) continue;
		apply_homography(p, iHa, m);
		apply_homography(q, iHb, m + 2);
		int ia = lrint(p[0]);
		int ja = lrint(p[1]);
		int ib = lrint(q[0]);
		int jb = lrint(q[1]);
		if (!insideP(wa, ha, ia, ja)) continue;
		if (!insideP(wb, hb, ib, jb)) continue;
		float *va = a + pda * (wa * ja + ia);
		float *vb = b + pda * (wb * jb + ib);
		if (!isfinite(*va)) continue;
		if (!isfinite(*vb)) continue;
		if (pda == 1) // heights
			fprintf(fo, "%lf %lf %lf %lf %lf %lf\n",
				p[0], p[1], *va, q[0], q[1], *vb);
		else // xyz
			fprintf(fo, "%lf %lf %lf %lf %lf %lf\n",
				va[0], va[1], va[2], vb[0], vb[1], vb[2]);
	}

	free(a);
	free(b);
	xfclose(fo);
	xfclose(fi);
	return 0;
}
