#include "vvector.h"


// RANSAC INPUT EXAMPLES

// ************************************
// ************************************
// ************************************
// ************************************
//
// 1. straight lines
// 	datadim = 2 (coordinates of each point)
// 	modeldim = 3 (coefficients of straight line)
// 	nfit = 2 (number of points that determine a straight line)

// instance of "ransac_error_evaluation_function"
float distance_of_point_to_straight_line(float *line, float *point, void *usr)
{
	float n = hypot(line[0], line[1]);
	float a = line[0]/n;
	float b = line[1]/n;
	float c = line[2]/n;
	float x = point[0];
	float y = point[1];
	float e = a*x + b*y + c;
	return fabs(e);
}


// instance of "ransac_model_generating_function"
int straight_line_through_two_points(float *line, float *points, void *usr)
{
	float ax = points[0];
	float ay = points[1];
	float bx = points[2];
	float by = points[3];
	float n = hypot(bx - ax, by - ay);
	if (!n) return 0;
	float dx = -(by - ay)/n;
	float dy = (bx - ax)/n;
	line[0] = dx;
	line[1] = dy;
	line[2] = -(dx*ax + dy*ay);

	// assert that the line goes through the two points
	float e1 = distance_of_point_to_straight_line(line, points, NULL);
	float e2 = distance_of_point_to_straight_line(line, points+2, NULL);
	assert(hypot(e1, e2) < 0.001);
	return 1;
}

int find_straight_line_by_ransac(bool *out_mask, float line[3],
		float *points, int npoints,
		int ntrials, float max_err)
{
	return ransac(out_mask, line, points, 2, npoints, 3,
			distance_of_point_to_straight_line,
			straight_line_through_two_points,
			4, ntrials, 3, max_err, NULL, NULL);
}

// ************************************
// ************************************
// ************************************
// ************************************
//
// 2. affine maps between pairs of points
// 	datadim = 4 (coordinates of each pair of points)
// 	modeldim = 6 (coefficients of the affine map)
// 	nfit = 3 (number of pairs that determine an affine map)

// utility function: evaluate an affine map
static void affine_map(float *y, float *A, float *x)
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

// instance of "ransac_error_evaluation_function"
static float affine_match_error(float *aff, float *pair, void *usr)
{
	float p[2] = {pair[0], pair[1]};
	float q[2] = {pair[2], pair[3]};
	float ap[2];
	affine_map(ap, aff, p);
	float e = hypot(ap[0] - q[0], ap[1] - q[1]);
	return e;
}


// naive linear algebra
static void apply_matrix6(double y[6], double A[6][6], double x[6])
{
	for (int i = 0; i < 6; i++)
	{
		double r = 0;
		for (int j = 0; j < 6; j++)
			r += A[i][j] * x[j];
		y[i] = r;
	}
}

// FIND AFFINITY
//
// Assuming that there is an affine tranformation
//
// A*x = x'
//
// or, in coordinates,
//
// /p q a\   /x\   /x'|
// |r s b| . |y| = |x'|
// \0 0 1/   \1/   \1 /
//
// the following function finds the matrix A from a set of three point pairs.
// It returns a condition number.  If the condition number is too close to
// zero, the result may be bad.
//
// A[] = {p, q, a, r, s, b}
//
static
double find_affinity(double *A, double *x, double *y, double *xp, double *yp)
{
	double caca = y[1]*x[0] + y[2]*x[1] + x[2]*y[0]
					- x[1]*y[0] - x[2]*y[1] - y[2]*x[0];
	double invB[6][6] = {
		         {y[1]-y[2], y[2]-y[0], y[0]-y[1], 0, 0, 0},
		         {x[2]-x[1], x[0]-x[2], x[1]-x[0], 0, 0, 0},
		         {x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2],
						x[0]*y[1]-x[1]*y[0], 0, 0, 0},
		{0, 0, 0, y[1]-y[2], y[2]-y[0], y[0]-y[1]},
		{0, 0, 0, x[2]-x[1], x[0]-x[2], x[1]-x[0]},
		{0, 0, 0, x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2],
						y[1]*x[0]-x[1]*y[0]}
	};
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			invB[j][i] /= caca;
	double Xp[6] = {xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]};
	apply_matrix6(A, invB, Xp);
	return caca;
}

// instance of "ransac_model_generating_function"
int affine_map_from_three_pairs(float *aff, float *pairs, void *usr)
{
	// call the function "find_affinity" from elsewhere
	double x[3] = {pairs[0], pairs[4], pairs[8]};
	double y[3] = {pairs[1], pairs[5], pairs[9]};
	double xp[3] = {pairs[2], pairs[6], pairs[10]};
	double yp[3] = {pairs[3], pairs[7], pairs[11]};
	double A[6];
	double r = find_affinity(A, x, y, xp, yp);
	for (int i = 0; i < 6; i++)
		aff[i] = A[i];
	return fabs(r) > 0.001;
}

// instance of "ransac_model_accepting_function"
bool affine_map_is_reasonable(float *aff, void *usr)
{
	float a = aff[0];
	float b = aff[1];
	float c = aff[3];
	float d = aff[4];
	float det = a*d - b*c;
	if (det < 0) return false;
	if (fabs(det) > 100) return false;
	if (fabs(det) < 0.01) return false;
	double n = a*a + b*b + c*c + d*d;
	if (n > 10) return false;
	if (n < 0.1) return false;
	return true;
}

int find_affinity_by_ransac(bool *out_mask, float affinity[6],
		float *pairs, int npairs,
		int ntrials, float max_err)
{
	return ransac(out_mask, affinity, pairs, 4, npairs, 6,
			affine_match_error,
			affine_map_from_three_pairs,
			3, ntrials, 4, max_err,
			affine_map_is_reasonable, NULL);
}

// ************************************
// ************************************
// ************************************
// ************************************
//
// 3. projective maps between pairs of points
// 	datadim = 4 (coordinates of each pair of points)
// 	modeldim = 9 (coefficients of the homographic map)
// 	nfit = 4 (number of pairs that determine an homographic map)

// ************************************
// ************************************
// ************************************
// ************************************
//
//
// TODO: homograpy, fundamental matrix

// utility function homographic map

// instance of "ransac_error_evaluation_function"
static float homographic_match_error(float *hom, float *pair, void *usr)
{
	double p[2] = {pair[0], pair[1]};
	double q[2] = {pair[2], pair[3]};
	const double H[3][3] = {{hom[0], hom[1], hom[2]},
		{hom[3], hom[4], hom[5]},
		{hom[6], hom[7], hom[8]}};
	double Hp[2];
	homography_transform(p, H, Hp);
	double r = hypot(Hp[0] - q[0], Hp[1] - q[1]);
	return r;
}

// instance of "ransac_model_generating_function"
int homography_from_four(float *hom, float *pairs, void *usr)
{
	double a[2] = {pairs[0], pairs[1]};
	double x[2] = {pairs[2], pairs[3]};
	double b[2] = {pairs[4], pairs[5]};
	double y[2] = {pairs[6], pairs[7]};
	double c[2] = {pairs[8], pairs[9]};
	double z[2] = {pairs[10], pairs[11]};
	double d[2] = {pairs[12], pairs[13]};
	double w[2] = {pairs[14], pairs[15]};
	double R[3][3], *RR=R[0];
	homography_from_4corresp(a, b, c, d, x, y, z, w, R);
	int r = 1;
	for (int i = 0; i < 9; i++)
		if (isfinite(RR[i]))
			hom[i] = RR[i];
		else
			r = 0;
	return r;
}


// ************************************
// ************************************
// ************************************
// ************************************
//
// 4. fundamental matrices
// 	datadim = 4 (coordinates of each pair)
// 	modeldim = 9 (fundamental matrix)
// 	nfit = 7 (seven-point algorithm)


static float fnorm(float *a, int n)
{
	if (n == 1)
		return abs(*a);
	if (n == 2)
		return hypot(a[0], a[1]);
	else
		return hypot(*a, fnorm(a+1, n-1));
}

#include "moistiv_epipolar.c"
// instance of "ransac_model_generating_function"
int seven_point_algorithm(float *fm, float *p, void *usr)
{
	int k[7] = {0, 1, 2, 3, 4, 5, 6};
	float m1[14] = {p[0],p[1], p[4],p[5], p[8],p[9],   p[12],p[13],
		          p[16],p[17], p[20],p[21], p[24],p[25] };
	float m2[14] = {p[2],p[3], p[6],p[7], p[10],p[11], p[14],p[15],
		          p[18],p[19], p[22],p[23], p[26],p[27] };
	float z[3], F1[4][4]={0}, F2[4][4]={0};
	// this is so braindead it's not funny
	int r = moistiv_epipolar(m1, m2, k, z, F1, F2);
	//MAT_PRINT_4X4(F1);
	//MAT_PRINT_4X4(F2);
	int ridx = 0;
	if (r == 3) ridx = random_index(0, 3);
	//fprintf(stderr, "z = %g\n", z[ridx]);
	int cx = 0;
	for (int j = 1; j <= 3; j++)
	for (int i = 1; i <= 3; i++)
		fm[cx++] = F1[i][j] + z[ridx]*F2[i][j];
	//fprintf(stderr, "\tspa = %g %g %g   %g %g %g   %g %g %g\n",
	//		fm[0], fm[1], fm[2],
	//		fm[3], fm[4], fm[5],
	//		fm[6], fm[7], fm[8]);
	float nfm = fnorm(fm, 9);
	for (int i = 0; i < 9; i++)
		fm[i] /= nfm;
	return 1;
}

// instance of "ransac_error_evaluation_function"
static float epipolar_algebraic_error(float *fm, float *pair, void *usr)
{
	float p[3] = {pair[0], pair[1], 1};
	float q[3] = {pair[2], pair[3], 1};
	float fp[3] = {fm[0]*p[0] + fm[1]*p[1] + fm[2],
		       fm[3]*p[0] + fm[4]*p[1] + fm[5],
		       fm[6]*p[0] + fm[7]*p[1] + fm[8]};
	float qfp = q[0]*fp[0] + q[1]*fp[1] + q[2]*fp[2];
	return fabs(qfp);
}

// instance of "ransac_error_evaluation_function"
static float epipolar_euclidean_error(float *fm, float *pair, void *usr)
{
	float p[3] = {pair[0], pair[1], 1};
	float q[3] = {pair[2], pair[3], 1};
	float pf[3] = {fm[0]*p[0] + fm[3]*p[1] + fm[6],
		       fm[1]*p[0] + fm[4]*p[1] + fm[7],
		       fm[2]*p[0] + fm[5]*p[1] + fm[8]};
	float npf = hypot(pf[0], pf[1]);
	pf[0] /= npf; pf[1] /= npf; pf[2] /= npf;
	float pfq = pf[0]*q[0] + pf[1]*q[1] + pf[2]*q[2];
	return fabs(pfq);
}

// instance of "ransac_error_evaluation_function"
static float epipolar_symmetric_euclidean_error(float *fm, float *pair, void *u)
{
	fail("must transpose F!");
	float rpair[4] = {pair[2], pair[3], pair[0], pair[1]};
	float a = epipolar_euclidean_error(fm, pair, u);
	float b = epipolar_euclidean_error(fm, rpair, u);
	return hypot(a, b);
}

// instance of "ransac_error_evaluation_function"
static float epipolar_error(float *fm, float *pair, void *usr)
{
	ransac_error_evaluation_function *f;
	f = epipolar_euclidean_error;
	return f(fm, pair, usr);
}


static void mprod33(float *ab_out, float *a_in, float *b_in)
{
	float (*about)[3] = (void*)ab_out;
	float (*a)[3] = (void*)a_in;
	float (*b)[3] = (void*)b_in;
	float ab[3][3] = {0};
	for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 3; k++)
			ab[i][j] += a[i][k] * b[k][j];
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
		about[i][j] = ab[i][j];
}

int find_fundamental_matrix_by_ransac(bool *out_mask, float out_fm[9],
		float *pairs, int npairs,
		int ntrials, float max_err)
{
	// compute input statistics
	double pmean[4] = {0, 0, 0, 0}, pdev[4] = {0, 0, 0, 0};
	for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 4; j++)
		pmean[j] += (pairs[4*i+j])/npairs;
	for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 4; j++)
		pdev[j] += (fabs(pairs[4*i+j]-pmean[j]))/npairs;
	float PDev = (pdev[0]+pdev[1]+pdev[2]+pdev[3])/4;
	float pDev[4] = {PDev, PDev, PDev, PDev};

	fprintf(stderr, "pmean = %g %g %g %g\n",
			pmean[0], pmean[1], pmean[2], pmean[3]);
	fprintf(stderr, "pdev = %g %g %g %g {%g}\n",
			pdev[0], pdev[1], pdev[2], pdev[3], PDev);

	// normalize input
	float *pairsn = xmalloc(4*npairs*sizeof*pairsn);
	for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 4; j++)
		pairsn[4*i+j] = (pairs[4*i+j]-pmean[j])/PDev;

	for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 4; j++)
		fprintf(stderr, "%g %g %g %g\n",
				pairsn[4*i+0], pairsn[4*i+1],
				pairsn[4*i+2], pairsn[4*i+3]);
	float max_err_n = max_err / PDev;

	// set up algorithm context
	int datadim = 4;
	int modeldim = 9;
	int nfit = 7;
	ransac_error_evaluation_function *f_err = epipolar_error;
	ransac_model_generating_function *f_gen = seven_point_algorithm;
	ransac_model_accepting_function  *f_acc = NULL;

	// run algorithm on normalized data
	float nfm[9];
	int n_inliers = ransac(out_mask, nfm,
			pairsn, datadim, npairs, modeldim,
			f_err, f_gen, nfit, ntrials, nfit+1, max_err_n,
			f_acc, NULL);

	// un-normalize result
	float Atrans[9] = {1/pDev[0], 0, 0,
		           0, 1/pDev[1], 0,
			   -pmean[0]/pDev[0], -pmean[1]/pDev[1], 1};
	float B[9] = {1/pDev[2], 0, -pmean[2]/pDev[2],
		      0, 1/pDev[3], -pmean[3]/pDev[3],
		      0, 0, 1};
	float tmp[9];
	mprod33(tmp, nfm, B);
	mprod33(out_fm, Atrans, tmp);

	// print matrices
	fprintf(stderr, "Atrans= [%g %g %g ; %g %g %g ; %g %g %g]\n",
			Atrans[0],Atrans[1],Atrans[2],Atrans[3],Atrans[4],Atrans[5],Atrans[6],Atrans[7],Atrans[8]);
	fprintf(stderr, "nfm= [%g %g %g ; %g %g %g ; %g %g %g]\n",
			nfm[0],nfm[1],nfm[2],nfm[3],nfm[4],nfm[5],nfm[6],nfm[7],nfm[8]);
	fprintf(stderr, "B= [%g %g %g ; %g %g %g ; %g %g %g]\n",
			B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8]);
	fprintf(stderr, "out_fm= [%g %g %g ; %g %g %g ; %g %g %g]\n",
			out_fm[0],out_fm[1],out_fm[2],out_fm[3],out_fm[4],out_fm[5],out_fm[6],out_fm[7],out_fm[8]);

	// for each pair, display its normalized and un-normalized errors
	for (int i = 0; i < npairs; i++)
	{
		float en = f_err(nfm, pairsn+4*i, NULL);
		float en2 = epipolar_algebraic_error(nfm, pairsn+4*i, NULL);
		float eu = f_err(out_fm, pairs+4*i, NULL);
		fprintf(stderr, "en,eu = %g\t%g\t{%g}\t%g\n", en, eu,eu/en,en2);
	}


	free(pairsn);
	return n_inliers;
}
