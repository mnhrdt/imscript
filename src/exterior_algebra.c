// NOTE: for odd-dimensional spaces, the bases are chosen so that the Hodge
// linear map is defined by the identity matrix.

// general code to compute an exterior product defined by a table t
void exterior_product(double *xy, double *x, double *y, int (*t)[4],
		int n_xy, int n_t, int fold)
{
	// initialize the output with zeros
	for (int i = 0; i < n_xy; i++)
		xy[i] = 0;

	// accumulate
	for (int i = 0; i < n_t; i++)
	{
		int a = t[i][0]; // basis index for x
		int b = t[i][1]; // basis index for y
		int c = t[i][2]; // basis index for xy
		int k = t[i][3]; // coefficient
		xy[c] += k * x[a] * y[b];
		if (fold)
			xy[c] += fold * k * x[b] * y[a];
	}
}

// exterior algebra in R^3
// O^0 = < 1 >
// O^1 = < dx, dy, dz >
// O^2 = < dyz, dzx, dxy >
// O^3 = < dxyz >
enum { D3_x , D3_y , D3_z  };
enum { D3_yz, D3_zx, D3_xy };
enum { D3_xyz              };

// ^ : O^1 x O^1 -> O^2  (R^3 x R^3 -> R^3)
void exterior_product_3_3_3(double xy[3], double x[3], double y[3])
{
	int t[][4] = {
		{ D3_x, D3_y, D3_xy, 1},
		{ D3_x, D3_z, D3_zx, -1},
		{ D3_y, D3_z, D3_yz, 1}
	};
	exterior_product(xy, x, y, t, 3, 3, -1);
}

// ^ : O^1 x O^2 -> O^3  (R^3 x R^3 -> R)
void exterior_product_1_3_3(double xy[1], double x[3], double y[3])
{
	int t[][4] = {
		{ D3_x, D3_yz, D3_xyz, 1},
		{ D3_y, D3_zx, D3_xyz, 1},
		{ D3_z, D3_xy, D3_xyz, 1}
	}; // note: this is just the regular scalar product
	exterior_product(xy, x, y, t, 1, 3, 0);
}

// * : O^1 -> O^2 (using euclidean metric)
void exterior_star_3(double y[3], double x[3])
{
	// in the cosen order of the basis, this is just the identity
	y[0] = x[0];
	y[1] = x[1];
	y[2] = x[2];
}


// exterior algebra in R^4 (using minkowski metrix)
// O^0 = (1) < 1 >
// O^1 = (4) < dx, dy, dz, dt >
// O^2 = (6) < dyz, dzx, dxy, dxt, dyt, dzt >
// O^3 = (4) < dyzt, dzxt, dxyt, dxyz >
// O^4 = (1) < dxyzt >

enum { D4_x  , D4_y  , D4_z  , D4_t   }; // 4
enum { D4_yzt, D4_zxt, D4_xyt, D4_xyz }; // 4
enum { D4_yz, D4_zx, D4_xy, D4_xt, D4_yt, D4_zt }; // 6
enum { D4_xyzt };  // 1

// ^ : O^1 x O^1 -> O^2  (R^4 x R^4 -> R^6)
void exterior_product_6_4_4(double xy[6], double x[4], double y[4])
{
	int t[][4] = {
		{ D4_x, D4_y, D4_xy, 1},
		{ D4_x, D4_z, D4_zx, -1},
		{ D4_x, D4_t, D4_xt, 1},
		{ D4_y, D4_z, D4_yz, 1},
		{ D4_y, D4_t, D4_yt, 1},
		{ D4_z, D4_t, D4_zt, 1}
	};
	exterior_product(xy, x, y, t, 6, 6, -1);
}

// ^ : O^1 x O^2 -> O^3  (R^4 x R^6 -> R^4)
void exterior_product_4_4_6(double xy[4], double x[4], double y[6]);

// ^ : O^1 x O^3 -> O^4  (R^4 x R^4 -> R  )
void exterior_product_1_4_4(double xy[1], double x[4], double y[4]);

// ^ : O^2 x O^2 -> O^4  (R^6 x R^6 -> R  )
void exterior_product_1_6_6(double x[6], double y[6]);

// * : O^6 -> O^6
void exterior_star_4_2(double xy[1], double x[6], double y[6]);

// exterior algebra in R^5
// O^0 = (1)  < 1 >
// O^1 = (5)  < dx   , dy   , dp   , dq    , dt >
// O^2 = (10) < dxy, dxt, dyt, dpq, dpt, dqt, dxp, dxq, dyp, dyq >
// O^3 = (10) < ... >
// O^4 = (5)  < ... >
// O^5 = (1)  < dxypqt >
//

// The beauty of this algorithm is that *all* the complexity is hidden
// into variable names.  The values of these variables are arbitrary integers.
enum { D5_x   , D5_y   , D5_p   , D5_q   , D5_t    };
enum { D5_ypqt, D5_xqpt, D5_xyqt, D5_yxpt, D5_xypq };
enum { D5_xy ,D5_xp ,D5_xq ,D5_xt ,D5_yp ,D5_yq ,D5_yt ,D5_pq ,D5_pt ,D5_qt  };
enum { D5_pqt,D5_qyt,D5_ypt,D5_ypq,D5_xqt,D5_xpt,D5_xpq,D5_xyt,D5_xyq,D5_xyp };
enum { D5_xypqt };

// ^ : O^1 x O^1 -> O^2  (R^5  x R^5   -> R^10)
void exterior_product_10_5_5(double xy[10], double x[5], double y[5])
{
	int t[][4] = {
		{ D5_x, D5_y, D5_xy, 1 },
		{ D5_x, D5_p, D5_xp, 1 },
		{ D5_x, D5_q, D5_xq, 1 },
		{ D5_x, D5_t, D5_xt, 1 },
		{ D5_y, D5_p, D5_yp, 1 },
		{ D5_y, D5_q, D5_yq, 1 },
		{ D5_y, D5_t, D5_yt, 1 },
		{ D5_p, D5_q, D5_pq, 1 },
		{ D5_p, D5_t, D5_pt, 1 },
		{ D5_q, D5_t, D5_qt, 1 }
	};
	exterior_product(xy, x, y, t, 10, 10, -1);
}

// ^ : O^1 x O^2 -> O^3  (R^5  x R^10  -> R^10)
void exterior_product_10_5_10(double xy[10], double x[5], double y[10]);

// ^ : O^1 x O^3 -> O^4  (R^5  x R^10  -> R^5 )
void exterior_product_5_5_10(double xy[5], double x[5], double y[10]);

// ^ : O^2 x O^2 -> O^4  (R^10 x R^10  -> R^5 )
void exterior_product_5_10_10(double xy[5], double x[10], double y[10])
{
	// Tabley McTableface
	int t[][4] = {
		{ D5_xy, D5_pq, D5_xypq,  1},
		{ D5_xy, D5_pt, D5_yxpt, -1},
		{ D5_xy, D5_qt, D5_xyqt,  1},
		{ D5_xp, D5_yq, D5_xypq, -1},
		{ D5_xp, D5_yt, D5_yxpt,  1},
		{ D5_xp, D5_qt, D5_xqpt, -1},
		{ D5_xq, D5_yp, D5_xypq,  1},
		{ D5_xq, D5_yt, D5_xyqt, -1},
		{ D5_xq, D5_pt, D5_xqpt,  1},
		{ D5_xt, D5_yp, D5_yxpt, -1},
		{ D5_xt, D5_yq, D5_xyqt,  1},
		{ D5_xt, D5_pq, D5_xqpt, -1},
		{ D5_yp, D5_qt, D5_ypqt,  1},
		{ D5_yq, D5_pt, D5_ypqt, -1},
		{ D5_yt, D5_pq, D5_ypqt,  1},
	};
	exterior_product(xy, x, y, t, 5, 15, 1);
}

// ^ : O^1 x O^4 -> O^5  (R^5  x R^5   -> R   )
void exterior_product_1_5_5(double xy[1], double x[5], double y[5])
{
	int t[][4] = {
		{ D5_x, D5_ypqt, D5_xypqt, 1 },
		{ D5_y, D5_xqpt, D5_xypqt, 1 },
		{ D5_p, D5_xyqt, D5_xypqt, 1 },
		{ D5_q, D5_yxpt, D5_xypqt, 1 },
		{ D5_t, D5_xypq, D5_xypqt, 1 }
	};
	exterior_product(xy, x, y, t, 1, 5, 0);
}

// ^ : O^2 x O^3 -> O^5  (R^10 x R^10  -> R   )
void exterior_product_1_10_10(double xy[1], double x[10], double y[10]);

//
//
// end of functions
//
//

//#define EXTERIOR_ALGEBRA_MAIN
#ifdef EXTERIOR_ALGEBRA_MAIN
// silly consistency checks
#include <stdio.h>
#include "random.c"
void cross_product(double xy[3], double x[3], double y[3])
{
	xy[0] = x[1]*y[2] - x[2]*y[1];
	xy[1] = x[2]*y[0] - x[0]*y[2];
	xy[2] = x[0]*y[1] - x[1]*y[0];
}
double inner_product(double *x, double *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}
int main(void)
{
	// check that the 3-dimensional exterior products
	// are indeed the traditional vector and scalar products
	double x[3], y[3], z[3];
	for (int i = 0; i < 3; i++)
	{
		x[i] = 100*random_uniform();
		y[i] = 100*random_uniform();
	}
	printf("x = %g %g %g\n", x[0], x[1], x[2]);
	printf("y = %g %g %g\n", y[0], y[1], y[2]);
	printf("\n");
	exterior_product_3_3_3(z, x, y);
	printf("x = %g %g %g\n", x[0], x[1], x[2]);
	printf("y = %g %g %g\n", y[0], y[1], y[2]);
	printf("z = %g %g %g\n", z[0], z[1], z[2]);
	printf("\n");
	cross_product(z, x, y);
	printf("z = %g %g %g\n", z[0], z[1], z[2]);
	double a[1], b[1];
	exterior_product_1_3_3(a, x, z);
	exterior_product_1_3_3(b, z, y);
	printf("a b = %g %g\n", *a, *b);
	printf("\n");
	printf("\n");

	// check that the 5-dimensional products are consistent
	double X[5][5] = {
		{ 1, 0, 0, 0, 1 },
		{ 0, 1, 0, 0, 0 },
		{ 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 1, 0 },
		{ 0, 0, 0, 0, 0 }
	};
	double Z[2][10];
	for (int j = 0; j < 4; j++)
	for (int i = 0; i < 5; i++)
		//X[j][i] = j==i;
		X[j][i] = 100*random_uniform();
printf("X[0] = %g %g %g %g %g\n",X[0][0],X[0][1],X[0][2],X[0][3],X[0][4]);
printf("X[1] = %g %g %g %g %g\n",X[1][0],X[1][1],X[1][2],X[1][3],X[1][4]);
printf("X[2] = %g %g %g %g %g\n",X[2][0],X[2][1],X[2][2],X[2][3],X[2][4]);
printf("X[3] = %g %g %g %g %g\n",X[3][0],X[3][1],X[3][2],X[3][3],X[3][4]);
	printf("\n");
	exterior_product_10_5_5(Z[0], X[0], X[1]);
	exterior_product_10_5_5(Z[1], X[2], X[3]);
	printf("Z[0] = %g %g %g %g %g %g %g %g %g %g\n",
			Z[0][0], Z[0][1], Z[0][2], Z[0][3], Z[0][4],
			Z[0][5], Z[0][6], Z[0][7], Z[0][8], Z[0][9]
	      );
	printf("Z[1] = %g %g %g %g %g %g %g %g %g %g\n",
			Z[1][0], Z[1][1], Z[1][2], Z[1][3], Z[1][4],
			Z[1][5], Z[1][6], Z[1][7], Z[1][8], Z[1][9]
	      );
	printf("\n");
	exterior_product_5_10_10(X[4], Z[0], Z[1]);
printf("X[0] = %g %g %g %g %g\n",X[0][0],X[0][1],X[0][2],X[0][3],X[0][4]);
printf("X[1] = %g %g %g %g %g\n",X[1][0],X[1][1],X[1][2],X[1][3],X[1][4]);
printf("X[2] = %g %g %g %g %g\n",X[2][0],X[2][1],X[2][2],X[2][3],X[2][4]);
printf("X[3] = %g %g %g %g %g\n",X[3][0],X[3][1],X[3][2],X[3][3],X[3][4]);
printf("X[4] = %g %g %g %g %g\n",X[4][0],X[4][1],X[4][2],X[4][3],X[4][4]);
	printf("\n");
	double A[4];
	exterior_product_1_5_5(A+0, X[0], X[4]);
	exterior_product_1_5_5(A+1, X[1], X[4]);
	exterior_product_1_5_5(A+2, X[2], X[4]);
	exterior_product_1_5_5(A+3, X[3], X[4]);
	printf("A = %g %g %g %g\n", A[0], A[1], A[2], A[3]);
	printf("\n");
	A[0] = inner_product(X[0], X[4], 5);
	A[1] = inner_product(X[1], X[4], 5);
	A[2] = inner_product(X[2], X[4], 5);
	A[3] = inner_product(X[3], X[4], 5);
	printf("A = %g %g %g %g\n", A[0], A[1], A[2], A[3]);
	printf("\n");
}
#endif
