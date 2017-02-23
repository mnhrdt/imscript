#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)



#include "fragments.c"

static float planar_gaussian(float cx, float cy, float s, float x, float y)
{
	float r = hypot(x-cx, y-cy);
	return exp(-(r/s)*(r/s));
}

// sigma = average flow magnitude
// mu = spatial smoothness
static void fill_random_flow(float (**x)[2], int w, int h, float sgm, float mu)
{
	int nblobs = 1;
	float p[2*nblobs][3];
	FORI(2*nblobs) {
		p[i][0] = w*random_uniform();
		p[i][1] = h*random_uniform();
		p[i][2] = mu*sqrt(w*h)*(0.5+0.5*random_uniform());
	}
	FORL(nblobs) {
		float tsgm = sgm*(2*random_uniform()-1);
		FORJ(h) FORI(w) {
			float fx = 0;//sgm * random_normal();
			float fy = 0;//sgm * random_normal();
			int n = nblobs;
			fx += planar_gaussian(p[l][0], p[l][1], p[l][2], i, j);
			fy +=planar_gaussian(p[2*n-l-1][0],p[2*n-l-1][1],p[2*n-l-1][2],i,j);
			x[j][i][0] = tsgm*fx;
			x[j][i][1] = tsgm*fy;
		}
	}
	double m = 0;
	FORJ(h) FORI(w)
		m += hypot(x[j][i][0], x[j][i][1]);
	m /= w*h;
	FORJ(h) FORI(w) FORL(2)
		x[j][i][l] *= sgm/m;

}

static void fill_cidentity(float (**x)[2], int w, int h, float p[1])
{
	float r = *p;
	float cx = w/2.0;
	float cy = h/2.0;
	FORJ(h) FORI(w) {
		float fx = (i - cx)/r;
		float fy = (j - cy)/r;
		x[j][i][0] = fx;
		x[j][i][1] = fy;
	}
}

static void fill_traslation(float (**x)[2], int w, int h, float p[2])
{
	FORJ(h) FORI(w) FORL(2)
		x[j][i][l] = p[l];
}

static void invert_traslation(float it[2], float t[2])
{
	FORL(2) it[l] = -t[l];
}

static void fill_affinity(float (**x)[2], int w, int h, float p[6])
{
	FORJ(h) FORI(w) {
		float fx = p[0]*i + p[1]*j + p[2];
		float fy = p[3]*i + p[4]*j + p[5];
		x[j][i][0] = fx - i;
		x[j][i][1] = fy - j;
	}
}

static void invert_affinity(float invA[6], float A[6])
{
	float a, b, c, d, p, q;
	a=A[0]; b=A[1]; p=A[2];
	c=A[3]; d=A[4]; q=A[5];
	float det = a*d - b*c;
	invA[0] = d;
	invA[1] = -b;
	invA[2] = b*q-d*p;
	invA[3] = -c;
	invA[4] = a;
	invA[5] = c*p-a*q;
	FORI(6) invA[i] /= det;
}



#include "vvector.h"

static void invert_homography(float invH[9], float H[9])
{
	float h[3][3] = { {H[0], H[1], H[2]},
			{H[3], H[4], H[5]},
			{H[6], H[7], H[8]}};
	float det;
	float ih[3][3];
	INVERT_3X3(ih, det, h);
	FORI(9) invH[i] = ih[0][i];
}

static void affine_mapf(float y[2], float A[6], float x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

static void projective_mapf(float y[2], float H[9], float x[2])
{
	float z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
}

static void fill_homography(float (**x)[2], int w, int h, float p[9])
{
	FORJ(h) FORI(w) {
		float fx = p[0]*i + p[1]*j + p[2];
		float fy = p[3]*i + p[4]*j + p[5];
		float fz = p[6]*i + p[7]*j + p[8];
		x[j][i][0] = fx/fz - i;
		x[j][i][1] = fy/fz - j;
	}
}

// (x,y) |-> (x+(a*R^2), y+(a*R^2))
// R^2 = x*x + y*y
// x' = x + x*(a*x*x + a*y*y)
// y' = y + y*(a*x*x + a*y*y)
// X' = X + a*X*X*X

//#include <complex.h>
static double solvecubicspecial(double a, double b)
{
	long double x;
	long double r;
	if (a < 0) {
		//double complex aa;
		//((double *)&aa)[0] = a;
		//((double *)&aa)[1] = 0;
		//double complex bb;
		//((double *)&bb)[0] = b;
		//((double *)&bb)[1] = 0;
		//printf("aa = %g + %g __I__\n", creal(aa), cimag(aa));
		//printf("bb = %g + %g __I__\n", creal(bb), cimag(bb));
		//double complex xx, rr;
		//xx = cpow(csqrt((27*aa*bb*bb+4)/aa)/(2*sqrt(27)*aa)+bb/(2*aa),1.0/3);
		//printf("xx = %g + %g __I__\n", creal(xx), cimag(xx));
		//rr = xx-1/(3*aa*xx);
		//printf("rr = %g + %g __I__\n", creal(rr), cimag(rr));
		//return creal(rr);

		//error("not yet inverted");
		// the case a<0 corresponds to a cubic having three real roots,
		// and Cardano's formula involves complex numbers.  To avoid
		// complex numbers, the trigonometric representation of
		// solutions should be used instead.
// a*X*X*X + X - X' = 0
// X*X*X + (1/a)*X - X'/a = 0
// p = 1/a, q=b/a;
		long double p = 1/a;
		long double q = -b/a;
		int k = 1;
		long double cosarg = acos((3*q)/(2*p)*sqrt(-3/p))/3-k*2*M_PI/3;
		r = 2*sqrt(-p/3) * cos(cosarg);
		return r;
	} else {
		x = cbrt(sqrt((27*a*b*b+4)/a)/(2*sqrt(27)*a)+b/(2*a));
		r = x-1/(3*a*x);
		return r;
	}
//	double x = sqrt((27*a*b + 4)/a) / (2*sqrt(27)*a) + b/(2*a);
//	double y = pow(x, 1.0/3);
//	//double r = y - 1/(3*a*y);
//	double r = (3*a*y*y - 1)/(3*a*y);
//	return r;
}

// X' = X + a*X*X*X
// a*X*X*X + X - X' = 0
static double invertparabolicdistortion(double a, double xp)
{
	return solvecubicspecial(a, xp);
}

// X' = X + a*X*X*X
static double parabolicdistortion(double a, double x)
{
	return x + a*x*x*x;
}

static void fill_radialpol(float (**x)[2], int w, int h, float *p)
{
	int np = *p;
	if (np != 4 && np != 5)
		error("bad parametric distortion np = %d\n", np);
	float c[2] = {p[1], p[2]};
	float *coef = p + 3;
	float a = coef[0];
	FORJ(h) FORI(w) {
		float r = hypot(i - c[0], j - c[1]);
		float R = parabolicdistortion(a, r);
		if (r > 0) {
			x[j][i][0] = c[0] + (R/r)*(i - c[0]) - i;
			x[j][i][1] = c[1] + (R/r)*(j - c[1]) - j;
		} else {
			x[j][i][0] = 0;
			x[j][i][1] = 0;
		}
	}
}

static void viewflow(uint8_t (**y)[3], float (**x)[2], int w, int h, float m)
{
	FORJ(h) FORI(w) {
		float *v = x[j][i];
		double r = hypot(v[0], v[1]);
		r = r>m ? 1 : r/m;
		double a = atan2(v[1], -v[0]);
		a = (a+M_PI)*(180/M_PI);
		a = fmod(a, 360);
		double hsv[3], rgb[3];
		hsv[0] = a;
		hsv[1] = r;
		hsv[2] = r;
		void hsv_to_rgb_doubles(double*, double*);
		hsv_to_rgb_doubles(rgb, hsv);
		FORL(3)
			y[j][i][l] = 255*rgb[l];

	}
}

static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	//case 0: return interpolate_nearest(a, b, c, d, x, y);
	//case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	default: error("caca de vaca");
	}
	return -1;
}


#include "getpixel.c"

static void general_interpolate(float *result,
		int pd, float (**x)[pd], int w, int h, float p, float q,
		int m) // method
{
	//if (p < 0 || q < 0 || p+1 >= w || q+1 >= h) {
	//	FORL(pd) result[l] = 0;
	//} else {
		int ip = floor(p);
		int iq = floor(q);
		FORL(pd) {
			float a = getsample_0(x[0][0], w, h, pd, ip  , iq  , l);
			float b = getsample_0(x[0][0], w, h, pd, ip  , iq+1, l);
			float c = getsample_0(x[0][0], w, h, pd, ip+1, iq  , l);
			float d = getsample_0(x[0][0], w, h, pd, ip+1, iq+1, l);
			//float a = x[iq][ip][l];
			//float b = x[iq+1][ip][l];
			//float c = x[iq][ip+1][l];
			//float d = x[iq+1][ip+1][l];
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			//fprintf(stderr, "p%g q%g ip%d iq%d a%g b%g c%g d%g l%d v%g\n", p, q, ip, iq, a, b, c, d, l, v);
			result[l] = v;
		}
	//}
}


static void apply_affinity(int pd, float (**y)[pd], float (**x)[pd],
		int w, int h, float A[6])
{
	float invA[6]; invert_affinity(invA, A);
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		affine_mapf(q, invA, p);
		float val[pd];
		general_interpolate(val, pd, x, w, h, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

static void apply_homography(int pd, float (**y)[pd], float (**x)[pd],
		int w, int h, float H[9])
{
	float invH[9]; invert_homography(invH, H);
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		projective_mapf(q, invH, p);
		float val[pd];
		general_interpolate(val, pd, x, w, h, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

static void apply_radialpol(int pd, float (**y)[pd], float (**x)[pd],
		int w, int h, float *param)
{
	int np = *param;
	if (np != 4)
		error("bad parametric distortion np = %d\n", np);
	float c[2] = {param[1], param[2]};
	float a = param[3];
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		float r = hypot(p[0]-c[0], p[1]-c[1]);
		float ir = invertparabolicdistortion(a, r);
		// q = c + ir*(p-c)
		q[0] = c[0] + (ir/r) * (p[0] - c[0]);
		q[1] = c[1] + (ir/r) * (p[1] - c[1]);
		float val[pd];
		general_interpolate(val, pd, x, w, h, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

#include "cmphomod.c"

static double produce_homography(double H[9], int w, int h,
		char *homtype, double *v)
{
	if (0 == strcmp(homtype, "hom")) { // actual parameters
		FORI(9) H[i] = v[i];
	} else if (0 == strcmp(homtype, "hom4p")) {
		// absolute displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + v[0], 0 + v[1]},
			{w + v[2], 0 + v[3]},
			{0 + v[4], h + v[5]},
			{w + v[6], h + v[7]}
		};
		double R[3][3];
		homography_from_4corresp(
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3], R);
		FORI(9) H[i] = R[0][i];
	} else if (0 == strcmp(homtype, "hom4pr")) {
		// absolute displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + w*v[0], 0 + h*v[1]},
			{w + w*v[2], 0 + h*v[3]},
			{0 + w*v[4], h + h*v[5]},
			{w + w*v[6], h + h*v[7]}
		};
		double R[3][3];
		homography_from_4corresp(
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3], R);
		FORI(9) H[i] = R[0][i];
	} else if (0 == strcmp(homtype, "hom4prc")) {
		// percentual relative displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + w*v[0]/100, 0 + h*v[1]/100},
			{w + w*v[2]/100, 0 + h*v[3]/100},
			{0 + w*v[4]/100, h + h*v[5]/100},
			{w + w*v[6]/100, h + h*v[7]/100}
		};
		double R[3][3];
		homography_from_4corresp(
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3], R);
		FORI(9) H[i] = R[0][i];
	} else if (0 == strcmp(homtype, "hom16")) {
		// absolute coordinates of 4 point pairs
		double a[4][2] = {
			{v[0],v[1]},{v[2],v[3]},{v[4],v[5]},{v[6],v[7]}
		};
		double b[4][2] = {
			{v[8],v[9]},{v[10],v[11]},{v[12],v[13]},{v[14],v[15]}
		};
		double R[3][3];
		homography_from_4corresp( a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3], R);
		FORI(9) H[i] = R[0][i];
	} else if (0 == strcmp(homtype, "hom16r")) {
		// relative coordinates of 4 point pairs
		double a[4][2] = {
			{w*v[0],h*v[1]},
			{w*v[2],h*v[3]},
			{w*v[4],h*v[5]},
			{w*v[6],h*v[7]}
		};
		double b[4][2] = {
			{w*v[8], h*v[9]},
			{w*v[10],h*v[11]},
			{w*v[12],h*v[13]},
			{w*v[14],h*v[15]}
		};
		double R[3][3];
		homography_from_4corresp( a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3], R);
		FORI(9) H[i] = R[0][i];
	} else if (0 == strcmp(homtype, "hom16rc")) {
		// percentual relative coordinates of 4 point pairs
		double a[4][2] = {
			{w*v[0]/100,h*v[1]/100},
			{w*v[2]/100,h*v[3]/100},
			{w*v[4]/100,h*v[5]/100},
			{w*v[6]/100,h*v[7]/100}
		};
		double b[4][2] = {
			{w*v[8] /100,h*v[9] /100},
			{w*v[10]/100,h*v[11]/100},
			{w*v[12]/100,h*v[13]/100},
			{w*v[14]/100,h*v[15]/100}
		};
		double R[3][3];
		homography_from_4corresp( a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3], R);
		FORI(9) H[i] = R[0][i];
	} else error("unrecognized homography type \"%s\"", homtype);
	return 0;
}


// cid ZOOM
// cidb DISP
// tr DX DY
// aff A B P C D Q
// aff3p X1 Y1 X2 Y2 X3 Y3 X1' Y1' X2' Y2' X3' Y3'
// hom V1 ... V9
// hom4p X1 Y1 X2 Y2 X3 Y3 X4 Y4 X1' Y1' X2' Y2' X3' Y3' X4' Y4'
// hom4pr X1 Y1 X2 Y2 X3 Y3 X4 Y4 X1' Y1' X2' Y2' X3' Y3' X4' Y4'
// hom4prc X1 Y1 X2 Y2 X3 Y3 X4 Y4 X1' Y1' X2' Y2' X3' Y3' X4' Y4'
// hom16 p1 ... p16
// hom16r p1 ... p16
// hom16rc p1 ... p16
// hombord DX1 DY1 DX2 DY2 DX3 DY3 DX4 DY4
// smooth SIGMA MU
// cradial2 R
// cradial4 A B
// radial2 CX CY R
int main_synflow(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s {tr|aff|hom|smooth...} \"params\""
				//         0        1                   2
				" in out flow\n", *v);
				//3  4   5
		return EXIT_FAILURE;
	}

	int w, h, pd;
	void *data = iio_read_image_float_matrix_vec(v[3], &w, &h, &pd);
	float (**x)[pd] = data;

	float (**y)[pd] = NULL; y = matrix_build(w, h, sizeof**y);
	float (**f)[2] = matrix_build(w, h, sizeof**f);
	uint8_t (**vf)[3] = matrix_build(w, h, sizeof**vf);

	int maxparam = 40;
	float param[maxparam];
	int nparams = parse_floats(param, maxparam, v[2]);
	double dparam[maxparam]; FORI(nparams) dparam[i] = param[i];
	if (false) { ;
	} else if (0 == strcmp(v[1], "cradial2")) {
		if (nparams != 1) error("\"cradial2\" expects one parameter");
		float p[4] = {4, w/2.0, h/2.0, param[0]};
		fill_radialpol(f, w, h, p);
		apply_radialpol(pd, y, x, w, h, p);
		//double a = -0.01;
		//FORI(11) {
		//	double X = i/10.0;
		//	double xp = parabolicdistortion(a, X);
		//	double ixp = invertparabolicdistortion(a, xp);
		//	printf("x=%g, xp=%g, ixp=%g\n", X, xp, ixp);
		//}
		//return EXIT_SUCCESS;
	} else if (0 == strcmp(v[1], "cid")) {
		if (nparams != 1) error("\"cid\" expects one parameter");
		//fill_cidentity(f, w, h, param);
		float R = param[0];
		float tp[6] = {R, 0, (1-R)*w/2.0, 0, R, (1-R)*h/2.0};
		fill_affinity(f, w, h, tp);
		apply_affinity(pd, y, x, w, h, tp);
	} else if (0 == strcmp(v[1], "tr")) {
		if (nparams != 2) error("\"tr\" expects two parameters");
		//fill_traslation(f, w, h, param);
		float tp[6] = {1, 0, param[0], 0, 1, param[1]};
		fill_affinity(f, w, h, tp);
		apply_affinity(pd, y, x, w, h, tp);
	} else if (0 == strcmp(v[1], "aff")) {
		if (nparams != 6) error("\"aff\" expects six parameters");
		fill_affinity(f, w, h, param);
		apply_affinity(pd, y, x, w, h, param);
	} else if
		(v[1] == strstr(v[1], "hom")) {
		double Hd[9];
		float Hf[9];
		produce_homography(Hd, w, h, v[1], dparam);
		FORI(9) Hf[i] = Hd[i];
		fill_homography(f, w, h, Hf);
		apply_homography(pd, y, x, w, h, Hf);
//		(0 == strcmp(v[1], "hom")) {
//		if (nparams != 9) error("\"hom\" expects nine parameters");
//		fprintf(stderr, "applying homography %g %g %g %g %g %g %g %g %g"
//				" \n",
//				param[0], param[1], param[2],
//				param[3], param[4], param[5],
//				param[6], param[7], param[8]);
//		fill_homography(f, w, h, param);
//		apply_homography(pd, y, x, w, h, param);
//	} else if (0 == strcmp(v[1], "hom4pr")) {
//		if (nparams != 9) error("\"hom4pr\" expects eight parameters");
//		fprintf(stderr, "applying hom4pr %g %g %g %g %g %g %g %g"
//				" \n",
//				param[0], param[1], param[2],
//				param[3], param[4], param[5],
//				param[6], param[7]);
//		float H[9];
//		compute_h
//		fill_homography(f, w, h, H);
//		apply_homography(pd, y, x, w, h, param);
	} else if (0 == strcmp(v[1], "smooth")) {
		if (nparams != 2) error("\"smooth\" expects two parameters");
		fill_random_flow(f, w, h, param[0], param[1]);
		error("not yet implemented");
	} else
		error("unrecognized option \"%s\"\n", v[1]);


	//viewflow(vf, f, w, h, 1);
	//iio_write_image_uint8_vec("/tmp/vf.png", vf[0][0], w, h, 3);

	iio_write_image_float_vec(v[4], y[0][0], w, h, pd);
	iio_write_image_float_vec(v[5], f[0][0], w, h, 2);

	free(x);
	free(y);
	free(f);
	free(vf);

	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_synflow(c, v);
}
#endif//OMIT_MAIN
