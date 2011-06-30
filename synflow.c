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


static void general_interpolate(float *result,
		int pd, float (**x)[pd], int w, int h, float p, float q,
		int m) // method
{
	if (p < 0 || q < 0 || p+1 >= w || q+1 >= h) {
		FORL(pd) result[l] = 0;
	} else {
		int ip = floor(p);
		int iq = floor(q);
		FORL(pd) {
			//float a = getsample(x, w, h, pd, ip  , iq  , l);
			//float b = getsample(x, w, h, pd, ip  , iq+1, l);
			//float c = getsample(x, w, h, pd, ip+1, iq  , l);
			//float d = getsample(x, w, h, pd, ip+1, iq+1, l);
			float a = x[iq][ip][l];
			float b = x[iq+1][ip][l];
			float c = x[iq][ip+1][l];
			float d = x[iq+1][ip+1][l];
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			//fprintf(stderr, "p%g q%g ip%d iq%d a%g b%g c%g d%g l%d v%g\n", p, q, ip, iq, a, b, c, d, l, v);
			result[l] = v;
		}
	}
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

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}


// cid ZOOM
// cidb DISP
// tr DX DY
// aff A B P C D Q
// aff3p X1 Y1 X2 Y2 X3 Y3 X1' Y1' X2' Y2' X3' Y3'
// hom V1 ... V9
// hom4p X1 Y1 X2 Y2 X3 Y3 X4 Y4 X1' Y1' X2' Y2' X3' Y3' X4' Y4'
// hombord DX1 DY1 DX2 DY2 DX3 DY3 DX4 DY4
// smooth SIGMA MU
// parabolic R S
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

	int maxparam = 10;
	float param[maxparam];
	int nparams = parse_floats(param, maxparam, v[2]);
	if (0 == strcmp(v[1], "cid")) {
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
	} else if (0 == strcmp(v[1], "hom")) {
		if (nparams != 9) error("\"hom\" expects nine parameters");
		fprintf(stderr, "applying homography %g %g %g %g %g %g %g %g %g"
				" \n",
				param[0], param[1], param[2],
				param[3], param[4], param[5],
				param[6], param[7], param[8]);
		fill_homography(f, w, h, param);
		apply_homography(pd, y, x, w, h, param);
	} else if (0 == strcmp(v[1], "smooth")) {
		if (nparams != 2) error("\"smooth\" expects two parameters");
		fill_random_flow(f, w, h, param[0], param[1]);
		error("not yet implemented");
	} else
		error("unrecognized option \"%s\"\n", v[1]);


	//viewflow(vf, f, w, h, 1);
	//iio_save_image_uint8_vec("/tmp/vf.png", vf[0][0], w, h, 3);

	iio_save_image_float_vec(v[4], y[0][0], w, h, pd);
	iio_save_image_float_vec(v[5], f[0][0], w, h, 2);

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
