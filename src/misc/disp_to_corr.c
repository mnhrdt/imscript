#include <stdio.h>
#include <stdlib.h>


static double determinant_3x3(double a[9])
{
	return a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
}

static void invert_homography(double b[9], double a[9])
{
	double det = determinant_3x3(a);
	b[0] = (a[4]*a[8]-a[5]*a[7])/det;
	b[1] = (a[2]*a[7]-a[1]*a[8])/det;
	b[2] = (a[1]*a[5]-a[2]*a[4])/det;
	b[3] = (a[5]*a[6]-a[3]*a[8])/det;
	b[4] = (a[0]*a[8]-a[2]*a[6])/det;
	b[5] = (a[2]*a[3]-a[0]*a[5])/det;
	b[6] = (a[3]*a[7]-a[4]*a[6])/det;
	b[7] = (a[1]*a[6]-a[0]*a[7])/det;
	b[8] = (a[0]*a[4]-a[1]*a[3])/det;
}

static void apply_homography(double y[2], double H[9], double x[2])
{
	double z[3];
	z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
	z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
	z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = z[0]/z[2];
	y[1] = z[1]/z[2];
}

int disp_to_corr(float *out, double H1[9], double H2[9],
		float *disp, float *mask, int w, int h)
{
	double invH1[9]; invert_homography(invH1, H1);
	double invH2[9]; invert_homography(invH2, H2);

	int cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (mask && mask[j*w+i]) continue;
		double Aij[2], aij[2] = {i, j};
		double Bij[2], bij[2] = {i + disp[j*w+i], j};
		apply_homography(Aij, H1, aij);
		apply_homography(Bij, H2, bij);
		out[4*cx+0] = Aij[0];
		out[4*cx+1] = Aij[1];
		out[4*cx+2] = Bij[0];
		out[4*cx+3] = Bij[1];
		cx += 1;
	}

	return cx;
}

#include "iio.h"
#include "fail.c"
#include "xmalloc.c"

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}


int main(int c, char *v[])
{
	if (c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t"
				"%s \"H1\" \"H2\" disp [mask] > corr\n", *v);
		//                0   1      2    3     4
		return 1;
	}
	char *ascii_h1 = v[1];
	char *ascii_h2 = v[2];
	char *filename_disp = v[3];
	char *filename_mask = v[4];

	double H1[9], H2[9];
	if (9 != parse_doubles(H1, 9, ascii_h1)) fail("H1 should be 9 numbers");
	if (9 != parse_doubles(H2, 9, ascii_h2)) fail("H2 should be 9 numbers");

	// assumes the disparities are stored on a gray image
	int w[2], h[2];
	float *disp = iio_read_image_float(filename_disp, w, h);
	float *mask = NULL;
	if (filename_mask) {
		mask = iio_read_image_float(filename_disp, w+1, h+1);
		if (w[0] != w[1] || h[0] != h[1])
			fail("image and mask size mismatch (%d %d)!=(%d %d)",
					w[0], h[0], w[1], h[1]);
	}

	float *o = xmalloc(4 * *w * *h * sizeof*o);
	int n = disp_to_corr(o, H1, H2, disp, mask, *w, *h);

	for (int i = 0; i < n; i++)
		printf("%g %g %g %g\n", o[4*i], o[4*i+1], o[4*i+2], o[4*i+3]);

	free(o);
	free(disp);
	if (mask) free(mask);

	return 0;
}
