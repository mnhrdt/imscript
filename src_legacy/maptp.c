// map text points


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fail.c"
#include "xmalloc.c"
#include "parsenumbers.c"

static void apply_homography_2d(double y[2], double H[9], double x[2])
{
	double y0 = H[0]*x[0] + H[1]*x[1] + H[2];
	double y1 = H[3]*x[0] + H[4]*x[1] + H[5];
	double dd = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = y0/dd;
	y[1] = y1/dd;
}

static void apply_homography_3d(double y[3], double H[16], double x[3])
{
	double y0 = H[0]*x[0] + H[1]*x[1] + H[2]*x[2] + H[3];
	double y1 = H[4]*x[0] + H[5]*x[1] + H[6]*x[2] + H[7];
	double y2 = H[8]*x[0] + H[9]*x[1] + H[10]*x[2] + H[11];
	double dd = H[12]*x[0] + H[13]*x[1] + H[14]*x[2] + H[15];
	y[0] = y0/dd;
	y[1] = y1/dd;
	y[2] = y2/dd;
}

static int apply_hom2d(double *out, double *in, int nsamples,
		double *pars, int npars)
{
	if (npars != 9) fail("hom2d expects 9 parameters");
	int pdim = 2;
	int npoints = nsamples/pdim;
	for (int i = 0; i < npoints; i++)
		apply_homography_2d(out+pdim*i, pars, in+pdim*i);
	return npoints;
}

static int apply_hom3d(double *out, double *in, int nsamples,
		double *pars, int npars)
{
	if (npars != 16) fail("hom2d expects 16 parameters");
	int pdim = 3;
	int npoints = nsamples/pdim;
	for (int i = 0; i < npoints; i++)
		apply_homography_3d(out+pdim*i, pars, in+pdim*i);
	return npoints;
}

int main(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
				"%s {hom2d|hom3d} \"params\" <points.txt\n",*v);
		//                0  1              2
		return EXIT_FAILURE;
	}
	char *map_id = v[1];
	char *parstring = v[2];

	int (*f)(double*,double*,int,double*,int);
	if (false) ;
	else if (0 == strcmp(map_id, "hom2d")) f = apply_hom2d;
	else if (0 == strcmp(map_id, "hom3d")) f = apply_hom3d;
	else fail("unrecognized map \"%s\"", map_id);

	double params[0x100];
	int nparams = parse_doubles(params, 0x100, parstring);

	int nsamples;
	double *samples = read_ascii_doubles(stdin, &nsamples);

	double *outsamples = xmalloc(nsamples * sizeof*outsamples);

	int npoints = f(outsamples, samples, nsamples, params, nparams);
	int pdim = nsamples/npoints;
	for (int i = 0; i < nsamples; i++)
		printf("%lf%c", outsamples[i], (1+i)%pdim ? ' ':'\n');

	return EXIT_SUCCESS;
}
