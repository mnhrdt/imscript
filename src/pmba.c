// compute the reprojection error
// paired with "minimize.c", this is a poor man's bundle adjustment

#include <stdio.h>
#include <math.h>

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// compute the norm of a vector
static double vector_norm(double x[3])
{
	return hypot(hypot(x[0], x[1]), x[2]);
}

// compute the distance between two straight lines
// each line is represented by 6 numbers: a point and a direction vector
// if the straight lines are parallel, it returns NAN, instead of their distance
static double distance_between_two_straight_lines(double a[6], double b[6])
{
	double *u = a + 3;
	double *v = b + 3;
	double w[3]; vector_product(w, u, v);
	double nw = vector_norm(w);
	if (!isnormal(nw))
		return NAN;
	double amb[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
	return fabs(scalar_product(w, amb)/nw);
}

static double matrix_inversion(double ia[9], double a[9])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	ia[0] = ( a[4] * a[8] - a[5] * a[7] ) / det;
	ia[1] = ( a[2] * a[7] - a[1] * a[8] ) / det;
	ia[2] = ( a[1] * a[5] - a[2] * a[4] ) / det;
	ia[3] = ( a[5] * a[6] - a[3] * a[8] ) / det;
	ia[4] = ( a[0] * a[8] - a[2] * a[6] ) / det;
	ia[5] = ( a[2] * a[3] - a[0] * a[5] ) / det;
	ia[6] = ( a[3] * a[7] - a[4] * a[6] ) / det;
	ia[7] = ( a[1] * a[6] - a[0] * a[7] ) / det;
	ia[8] = ( a[0] * a[4] - a[1] * a[3] ) / det;
	return det;
}

static double matrix_times_vector(double y[3], double A[9], double x[3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	y[0] = A[0] * x[0]  +  A[1] * x[1]  +  A[2] * x[2];
	y[1] = A[3] * x[0]  +  A[4] * x[1]  +  A[5] * x[2];
	y[2] = A[6] * x[0]  +  A[7] * x[1]  +  A[8] * x[2];
}

static void matrix_product(double ab[9], double a[9], double b[9])
{
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		ab[3*i+j] = 0;
		for (int k = 0; k < 3; k++)
			ab[3*i+j] += a[3*i+k] * b[3*k+j];
	}
}

// find the straight line determined by a pixel on the given camera
static void from_pixel_to_line(double l[6], double P[12], double ij[2])
{
	// 1. decompose projection matrix P into rotation and translation parts
	//
	// 0 1 2 3
	// 4 5 6 7
	// 8 9 10 11
	double R[9] = {P[0], P[1], P[2], P[4], P[5], P[6], P[8], P[9], P[10]};
	double d[3] = {-P[3], -P[7], -P[11]};

	// 2. set-up output location
	double *center = l;
	double *direction = l + 3;

	// 3. compute camera center
	double iR[9];
	matrix_inversion(iR, R);
	matrix_times_vector(center, iR, d);

	// 4. compute line of sight of this pixel
	double ij1[3] = {ij[0], ij[1], 1};
	matrix_times_vector(direction, iR, ij1);

	//fprintf(stderr, "FPTL(%g %g) = %g %g %g  %g %g %g\n",
	//		ij[0], ij[1], l[0], l[1], l[2], l[3], l[4], l[5]);
}

static double compute_reprojection_error(int n, double *P, int nc, double *c)
{
	long double r = 0;
	for (int k = 0; k < nc; k++) // for each correspondence
	for (int j = 0; j < n; j++)
	for (int i = 0; i < j; i++)  // for each camera pair (i<j)
	{
		double *ck = c + 2*n*k; // numbers of the k-th correspondence
		double *ki = ck + 2*i;  // i-th point of the k-th correspondence
		double *kj = ck + 2*j;  // j-th point of the k-th correspondence
		double *Pi = P + 12*i;  // i-th projection matrix
		double *Pj = P + 12*j;  // j-th projection matrix
		if (!isfinite(ki[0] + ki[1] + kj[0] + kj[1]))
			continue;
		//fprintf(stderr, "ki=(%g %g) kj=(%g %g)\n", ki[0], ki[1], kj[0], kj[1]);
		double lin_ki[9], lin_kj[9];
		from_pixel_to_line(lin_ki, Pi, ki);
		from_pixel_to_line(lin_kj, Pj, kj);
		double e = distance_between_two_straight_lines(lin_ki, lin_kj);
		//fprintf(stderr, "e[%d][%d,%d] = %g\n", k, i, j, e);
		if (isfinite(e))
			r += e * e;
	}
	return r;
}

#include "xfopen.c"
#include "parsenumbers.c"

int main(int c, char *v[])
{
	// check and process input arguments
	if (c < 26) {
falla:
		fprintf(stderr, "usage:\n\t%s 2ncols.txt P1 ... P12n\n", *v);
		//                          0 1          2      c-1
		return c;
	}
	int n = (c - 2)/12;
	//fprintf(stderr, "n_views = %d\n", n);
	if (12*n + 2 != c) goto falla;
	char *filename_corr = v[1];
	double P[n][12];
	for (int i = 0; i < n; i++)
	for (int k = 0; k < 12; k++)
	{
		P[i][k] = atof(v[2 + 12*i + k]);
		//fprintf(stderr, "P[%d][%d] = %lf\n", i, k, P[i][k]);
	}

	// read point correspondances
	FILE *f = xfopen(filename_corr, "r");
	int ncorr;
	double *corr = read_ascii_doubles(f, &ncorr);
	if (ncorr % (2*n))
		return fprintf(stderr, "file \"%s\" must contain a multiple of "
			"%d numbers, but has %d\n", filename_corr, 2*n, ncorr);
	ncorr /= 2*n;
	//fprintf(stderr, "n_matches = %d\n", ncorr);
	xfclose(f);

	// compute and print reprojection error
	double err = compute_reprojection_error(n, P[0], ncorr, corr);
	printf("%lf\n", err);

	// cleanup and exit
	free(corr);
	return 0;
}
