// compute an "alternate" bundle adjustment cost (ie distance between rays)
// paired with "minimize.c", this is a poor man's bundle adjustment

#include <math.h>
#include <stdio.h>

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

static void matrix_times_vector(double y[3], double A[9], double x[3])
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

static void fill_rotation_matrix_from_angles(double out[9], double a[3])
{
    double c = cos(a[0]);
    double s = sin(a[0]);
    double Rx[9] = {1, 0, 0,
                    0, c, -s,
                    0, s, c};

    c = cos(a[1]);
    s = sin(a[1]);
    double Ry[9] = {c, 0, s,
                    0, 1, 0,
                    -s, 0, c};

    c = cos(a[2]);
    s = sin(a[2]);
    double Rz[9] = {c, -s, 0,
                    s, c, 0,
                    0, 0, 1};
    double tmp[9];
    matrix_product(tmp, Rx, Ry);
    matrix_product(out, tmp, Rz);
}

static void fill_rotation_matrix_from_quaternion(double *R, double r[4])
{
    double x = r[0];
    double y = r[1];
    double z = r[2];
    double w = r[3];

    *R++ = 1 - 2*y*y - 2*z*z;
    *R++ = 2*x*y - 2*z*w;
    *R++ = 2*x*z + 2*y*w;

    *R++ = 2*x*y + 2*z*w;
    *R++ = 1 - 2*x*x - 2*z*z;
    *R++ = 2*y*z - 2*x*w;

    *R++ = 2*x*z - 2*y*w;
    *R++ = 2*y*z + 2*x*w;
    *R++ = 1 - 2*x*x - 2*y*y;
}

static void compute_left_3x3_block_of_camera_matrix(double M[9], double r[10],
        double P[3])
{
    // r contains the reference parameters of the camera in that order:
    // f, px, py, cx, cy, cz, qx, qy, qz, qw
    // c is the center and q the orientation quaternion.
    double K[9] = {r[0], 0, r[1], 0, r[0], r[2], 0, 0, 1};
    double R[9];
    fill_rotation_matrix_from_quaternion(R, r+6);
    double Rp[9];
    fill_rotation_matrix_from_angles(Rp, P);
    double tmp[9];
    matrix_product(tmp, Rp, R);
    matrix_product(M, K, tmp);
}

// find the straight line determined by a pixel on the given camera
static void from_pixel_to_line(double l[6], double r[10], int nvar, double *P,
        double ij[2])
{
	// 1. set-up output location
	double *center = l;
	double *direction = l + 3;

	// 2. compute camera center
    center[0] = r[3];
    center[1] = r[4];
    center[2] = r[5];
    if (nvar == 6) {
        center[0] += P[0];
        center[1] += P[1];
        center[2] += P[2];
        P += 3;
    }

	// 3. compute line of sight of this pixel
    double M[9];
    compute_left_3x3_block_of_camera_matrix(M, r, P);
	double iM[9];
	matrix_inversion(iM, M);
	double ij1[3] = {ij[0], ij[1], 1};
	matrix_times_vector(direction, iM, ij1);

	//fprintf(stderr, "FPTL(%g %g) = %g %g %g  %g %g %g\n",
	//		ij[0], ij[1], l[0], l[1], l[2], l[3], l[4], l[5]);
}

// compute alternate bundle adjustment cost (ie distance between rays)
static double compute_aba_cost(int n, double *ref, int nvar, double *P, int nc,
        double *c)
{
	long double r = 0;
	for (int k = 0; k < nc; k++) // for each correspondence
	for (int j = 0; j < n; j++)
	for (int i = 0; i < j; i++)  // for each camera pair (i<j)
	{
		double *ck = c + 2*n*k; // numbers of the k-th correspondence
		double *ki = ck + 2*i;  // i-th point of the k-th correspondence
		double *kj = ck + 2*j;  // j-th point of the k-th correspondence
		double *ri = ref + 10*i;  // i-th camera reference parameters
		double *rj = ref + 10*j;  // j-th camera reference parameters
		double *Pi = P + nvar*i;  // i-th perturbation parameters
		double *Pj = P + nvar*j;  // j-th perturbation parameters
		if (!isfinite(ki[0] + ki[1] + kj[0] + kj[1]))
			continue;
		//fprintf(stderr, "ki=(%g %g) kj=(%g %g)\n", ki[0], ki[1], kj[0], kj[1]);
		double lin_ki[6], lin_kj[6];
		from_pixel_to_line(lin_ki, ri, nvar, Pi, ki);
		from_pixel_to_line(lin_kj, rj, nvar, Pj, kj);
		double e = distance_between_two_straight_lines(lin_ki, lin_kj);
		//fprintf(stderr, "e[%d][%d,%d] = %g\n", k, i, j, e);
		if (isfinite(e))
			r += e * e;
	}
	return r;
}


#include <stdio.h>
#include <stdbool.h>
#include "xfopen.c"
#include "parsenumbers.c"
#include "pickopt.c"

int main(int c, char *v[])
{
	// check and process input arguments
	bool do_normalize = pick_option(&c, &v, "n", 0);
	bool do_root = pick_option(&c, &v, "r", 0);
	double units = atof(pick_option(&c, &v, "u", "1"));
	bool no_translation = pick_option(&c, &v, "-no-translation", 0);
    int nvar = no_translation ? 3 : 6;

	if (c < (2 + 2 * (1 + nvar))) {
falla:
		//                         0    1         2
		fprintf(stderr, "usage:\n\t%s 2ncols.txt ref1.txt ... refn.txt"
                " [--no-translation] P1 ... P3n [ ... P6n]\n", *v);
        //fprintf(stderr, "c: %d, nvar: %d\n", c, nvar);
		return c;
	}

	int n = (c - 2) / (1 + nvar);
	if ((1 + nvar) * n + 2 != c) goto falla;
	//fprintf(stderr, "n_views = %d\n", n);
	char *filename_corr = v[1];
	double P[n][nvar];
	for (int i = 0; i < n; i++)
	for (int k = 0; k < nvar; k++)
	{
		P[i][k] = atof(v[2+n + nvar*i + k]);
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

    // read reference camera parameters
    double ref[n][10];
    for (int i = 0; i < n; i++)
        read_n_doubles_from_string(ref[i], v[2+i], 10);

	// compute and print reprojection error
    double err = compute_aba_cost(n, ref[0], nvar, P[0], ncorr, corr);
	if (do_normalize) err /= ncorr;
	if (do_root) err = sqrt(err);
	err *= units;
	printf("%lf\n", err);

	// cleanup and exit
	free(corr);
	return 0;
}
