// NOTE: for odd-dimensional spaces, the bases are chosen so that the Hodge
// linear map is defined by the identity matrix.

void exterior_product(double *xy, double *x, double *y, int (*t)[4],
		int n_xy, int n_t)
{
	// if n_t < 0, activate folding (duplicate the table with minus signs)
	int fold = 0;
	if (n_t < 0) {
		n_t = -n_t;
		fold = 1;
	}

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
			x[c] -= k * x[b] * y[a];
	}
}

// exterior algebra in R^3
// O^0 = < 1 >
// O^1 = < dx, dy, dz >
// O^2 = < dyz, dzx, dxy >
// O^3 = < dxyz >

enum basis_O1_R3 { D3_x , D3_y , D3_z  };
enum basis_O2_R3 { D3_yz, D3_zx, D3_xy };
enum basis_O3_R3 { D3_xyz              };

// ^ : O^1 x O^1 -> O^2  (R^3 x R^3 -> R^3)
void exterior_product_3_3_3(double xy[3], double x[3], double y[3])
{
	int t[][4] = {
		{ D3_x, D3_y, D3_xy, 1},
		{ D3_x, D3_z, D3_zx, -1},
		{ D3_y, D3_z, D3_yz, 1}
	};
	exterior_product(xy, x, y, t, 3, -3);
}

// ^ : O^1 x O^2 -> O^3  (R^3 x R^3 -> R)
void exterior_product_1_3_3(double xy[1], double x[3], double y[3])
{
	int t[][4] = {
		{ D3_x, D3_yz, D3_xyz, 1},
		{ D3_y, D3_zx, D3_xyz, 1},
		{ D3_z, D3_xy, D3_xyz, 1}
	};
	exterior_product(xy, x, y, t, 1, 3);
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

enum basis_O1_R4 { D4_x  , D4_y  , D4_z  , D4_t   }; // 4
enum basis_O3_R4 { D4_yzt, D4_zxt, D4_xyt, D4_xyz }; // 4
enum basis_O2_R4 { D4_yz, D4_zx, D4_xy, D4_xt, D4_yt, D4_zt }; // 6
enum basis_O4_R4 { D4_xyzt };  // 1

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
	exterior_product(xy, x, y, t, 6, -6);
}

// ^ : O^1 x O^2 -> O^3  (R^4 x R^6 -> R^4)
void exterior_product_4_4_6(double xy[4], double x[4], double y[6])
{
	//// NOT IMPLEMENTED
	//int t[][4] = {
	//	{ D4_x, D4_yz, D4_xyz, 1},
	//	{ D4_x, D4_zt, D4_xzt, 1},
	//	{ D4_y, D4_zt, D4_yzt, 1},
	//	{ D4_y, D4_xt, D4_yxt, 1},
	//	{ D4_z, D4_xt, D4_yxt, 1},
	//	{ D4_z, D4_xt, D4_yxt, 1},
	//};
	//exterior_product(xy, x, y, t, 4, 12);
}

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
//              dpqt dpqy dxpq dxyt dqxy dxyp dyqt dypt dxqt dxpt
// O^3 = (10) < ... >
// O^4 = (5)  < dypqt, dpxqt, dxyqt, dxytp, dxypq >
// O^5 = (1)  < dxypqt >
//

enum { D5_x, D5_y, D5_p, D5_q, D5_t };
enum { D5_xy ,D5_xp ,D5_xq ,D5_xt ,D5_yp ,D5_yq ,D5_yt ,D5_pq ,D5_pt ,D5_qt  };
enum { D5_pqt,D5_qyt,D5_ypt,D5_ypq,D5_xyt,D5_xpt,D5_xpq,D5_xyt,D5_xyq,D5_xyp };

#define D5_x D5_ypqt
#define D5_y D5_xpqt
#define D5_p D5_xyqt
#define D5_q D5_xypt
#define D5_t D5_xyqp
#define D5_xyqpt 0

//#define D5_x 0
//#define D5_y 1
//#define D5_p 2
//#define D5_q 3
//#define D5_t 4
//#define D5_xy 0
//#define D5_xp 1
//#define D5_xq 2
//#define D5_xt 3
//#define D5_yp 4
//#define D5_yq 5
//#define D5_yt 6
//#define D5_pq 7
//#define D5_pt 8
//#define D5_qt 9
//#define D5_xy D5_pqt
//#define D5_xp D5_qyt
//#define D5_xq D5_ypt
//#define D5_xt D5_ypq
//#define D5_yp D5_xyt
//#define D5_yq D5_xpt
//#define D5_yt D5_xpq
//#define D5_pq D5_xyt
//#define D5_pt D5_xyq
//#define D5_qt D5_xyp
//#define D5_x D5_ypqt
//#define D5_y D5_xpqt
//#define D5_p D5_xyqt
//#define D5_q D5_xypt
//#define D5_t D5_xyqp
//#define D5_xyqpt 0

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
	}
	exterior_product(xy, x, y, t, 10);
}

// ^ : O^1 x O^2 -> O^3  (R^5  x R^10  -> R^10)
void exterior_product_10_5_10(double xy[10], double x[5], double y[10]);

// ^ : O^1 x O^3 -> O^4  (R^5  x R^10  -> R^5 )
void exterior_product_5_5_10(double xy[5], double x[5], double y[10]);

// ^ : O^1 x O^4 -> O^5  (R^5  x R^5   -> R   )
void exterior_product_5_10_10(double xy[5], double x[10], double y[10])
{
	exterior_product(xy, x, y, t, 10);
}

// ^ : O^2 x O^2 -> O^4  (R^10 x R^10  -> R^5 )
void exterior_product_1_5_5(double xy[1], double x[5], double y[5]);

// ^ : O^2 x O^3 -> O^5  (R^10 x R^10  -> R   )
void exterior_product_1_10_10(double xy[1], double x[10], double y[10]);


//void exterior_product_10_5_5(double xy[10], double x[5], double y[5])
//{
//	xy[0] = x[0] * y[1] - y[0] * x[1];
//	xy[1] = x[0] * y[2] - y[0] * x[2];
//}
#endif
