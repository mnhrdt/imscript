// rational polynomial coefficient stuff

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xfopen.c"



// rational polynomial coefficients (and those of the inverse model)
struct rpc {
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];
};

// set all the values of an rpc model to NAN
static void nan_rpc(struct rpc *p)
{
	size_t nd = sizeof*p/sizeof(double);
	double *t = (double*)p;
	for (int i = 0; i < nd; i++)
		t[i] = NAN;
}

// like strcmp, but finds a needle
static int strhas(char *haystack, char *needle)
{
	char *r = strstr(haystack, needle);
	return r ? 0 : 1;
}

// returns the first character not in the prefix, or NULL if no prefix
static char *prefix(char *s, char *p)
{
	int i = 0;
	while (s[i] && p[i] && s[i] == p[i])
		i++;
	if (!s[i] || p[i])
		return NULL;
	else
		return s + i;
}

// if  s is the form "%s_%d"%(s,p), return %d bounded on [0,19]
static int pix(char *s, char *p)
{
	char *idx = prefix(s, p);
	if (!idx) return -1;
	int r = atoi(idx) - 1;
	if (r < 0) r = 0;
	if (r > 19) r = 19;
	return r;
}

// add a value to the rpc
static void add_tag_to_rpc(struct rpc *p, char *tag, double x)
{
	int t;
	if (false);
	else if (0 == strhas(tag, "SAMP_OFF"))          p->offset [0] = x;
	else if (0 == strhas(tag, "SAMP_SCALE"))        p->scale  [0] = x;
	else if (0 == strhas(tag, "LINE_OFF"))          p->offset [1] = x;
	else if (0 == strhas(tag, "LINE_SCALE"))        p->scale  [1] = x;
	else if (0 == strhas(tag, "HEIGHT_OFF"))        p->offset [2] = x;
	else if (0 == strhas(tag, "HEIGHT_SCALE"))      p->scale  [2] = x;
	else if (0 == strhas(tag, "LONG_OFF"))          p->ioffset[0] = x;
	else if (0 == strhas(tag, "LONG_SCALE"))        p->iscale [0] = x;
	else if (0 == strhas(tag, "LAT_OFF"))           p->ioffset[1] = x;
	else if (0 == strhas(tag, "LAT_SCALE"))         p->iscale [1] = x;
	else if (0 == strhas(tag, "HEIGHT_OFF"))        p->ioffset[2] = x;
	else if (0 == strhas(tag, "HEIGHT_SCALE"))      p->iscale [2] = x;
	else if (0 <= (t=pix(tag, "SAMP_NUM_COEFF_")))  p->numx   [t] = x;
	else if (0 <= (t=pix(tag, "SAMP_DEN_COEFF_")))  p->denx   [t] = x;
	else if (0 <= (t=pix(tag, "LINE_NUM_COEFF_")))  p->numy   [t] = x;
	else if (0 <= (t=pix(tag, "LINE_DEN_COEFF_")))  p->deny   [t] = x;
	else if (0 <= (t=pix(tag, "iSAMP_NUM_COEFF_"))) p->inumx  [t] = x;
	else if (0 <= (t=pix(tag, "iSAMP_DEN_COEFF_"))) p->idenx  [t] = x;
	else if (0 <= (t=pix(tag, "iLINE_NUM_COEFF_"))) p->inumy  [t] = x;
	else if (0 <= (t=pix(tag, "iLINE_DEN_COEFF_"))) p->ideny  [t] = x;
}

// process an input line
static double get_xml_tagged_number(char *tag, char *line)
{
	char buf[0x100], buf2[0x100];
	double x;
	int r = sscanf(line, " <%[^>]>%lf</%[^>]>", buf, &x, buf2);
	if (r == 3) {
		strcpy(tag, buf);
		return x;
	} else return NAN;
}

// read an XML file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename)
{
	nan_rpc(p);
	FILE *f = xfopen(filename, "r");
	int n = 0x100, o = 1;
	while (1) {
		char line[n], tag[n], *sl = fgets(line, n, f);;
		if (!sl) break;
		if (o && 0 == strhas(line, "<Inverse_Model>")) o = 0;
		tag[0] = 'i';
		double x = get_xml_tagged_number(tag+1, line);
		if (isfinite(x)) {
			//fprintf(stderr, "%s [%d]: %g\n", tag+o, o, x);
			add_tag_to_rpc(p, tag+o, x);
		}
	}
	xfclose(f);
	p->ioffset[2] = p->offset[2];
	p->iscale[2] = p->scale[2];
}

#define FORI(n) for (int i = 0; i < (n); i++)

void print_rpc(FILE *f, struct rpc *p, char *n)
{
	FORI(20) fprintf(f, "rpc(\"%s\") numx[%d] = %g\n", n, i, p->numx[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") denx[%d] = %g\n", n, i, p->denx[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") numy[%d] = %g\n", n, i, p->numy[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") deny[%d] = %g\n", n, i, p->deny[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") inumx[%d] = %g\n", n, i, p->inumx[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") idenx[%d] = %g\n", n, i, p->idenx[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") inumy[%d] = %g\n", n, i, p->inumy[i]);
	FORI(20) fprintf(f, "rpc(\"%s\") ideny[%d] = %g\n", n, i, p->ideny[i]);
	FORI(3)fprintf(stderr,"rpc(\"%s\") scale[%d] = %g\n",n,i,p->scale[i]);
	FORI(3)fprintf(stderr,"rpc(\"%s\") offset[%d] = %g\n",n,i,p->offset[i]);
	FORI(3)fprintf(stderr,"rpc(\"%s\") iscale[%d] = %g\n",n,i,p->iscale[i]);
	FORI(3)fprintf(stderr,"rpc(\"%s\") ioffset[%d]=%g\n",n,i,p->ioffset[i]);

}

// evaluate a polynomial of degree 3
double eval_pol20(double c[20], double x, double y, double z)
{
	//double col = x;
	//double lig = y;
	//double alt = z;
	//double m[20] = {1, lig, col, alt, lig*col,
	//	lig*alt, col*alt, lig*lig, col*col, alt*alt,
	//	col*lig*alt, lig*lig*lig, lig*col*col, lig*alt*alt, lig*lig*col,
	//	col*col*col, col*alt*alt, lig*lig*alt, col*col*alt, alt*alt*alt
	//};
	double m[20] = {1, y, x, z, x*y,
		y*z, x*z, y*y, x*x, z*z,
		x*y*z, y*y*y, y*x*x, y*z*z, y*y*x,
		x*x*x, x*z*z, y*y*z, x*x*z, z*z*z};
	double r = 0;
	for (int i = 0; i < 20; i++)
		r += c[i]*m[i];
	return r;
}

//// evaluate a quotient of 2 polynomials
//double eval_rpc40(double n[20], double d[20], double x, double y, double z)
//{
//}

// evaluate the normalized direct rpc model
static void eval_nrpc(double *result,
		struct rpc *p, double x, double y, double z)
{
	double numx = eval_pol20(p->numx, x, y, z);
	double denx = eval_pol20(p->denx, x, y, z);
	double numy = eval_pol20(p->numy, x, y, z);
	double deny = eval_pol20(p->deny, x, y, z);
	result[0] = numx/denx;
	result[1] = numy/deny;
	//fprintf(stderr, "\t\tnrpc{%p}(%g %g %g)=>(%g %g)\n", p, x, y, z, result[0], result[1]);

}

// evaluate the normalized inverse rpc model
static void eval_nrpci(double *result,
		struct rpc *p, double x, double y, double z)
{
	double numx = eval_pol20(p->inumx, y, x, z);
	double denx = eval_pol20(p->idenx, y, x, z);
	double numy = eval_pol20(p->inumy, y, x, z);
	double deny = eval_pol20(p->ideny, y, x, z);
	result[1] = numx/denx;
	result[0] = numy/deny;
	//fprintf(stderr, "\t\tnrpci{%p}(%g %g %g)=>",p,x,y,z);
	//fprintf(stderr, "(%g %g)\n", result[0], result[1]);
}

// evaluate the direct rpc model
static void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z)
{
	double nx = (x - p->offset[0])/p->scale[0];
	double ny = (y - p->offset[1])/p->scale[1];
	double nz = (z - p->offset[2])/p->scale[2];
	double tmp[2];
	eval_nrpc(tmp, p, nx, ny, nz);
	result[0] = tmp[0] * p->iscale[0] + p->ioffset[0];
	result[1] = tmp[1] * p->iscale[1] + p->ioffset[1];
}

// evaluate the inverse rpc model
static void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z)
{
	double nx = (x - p->ioffset[0])/p->iscale[0];
	double ny = (y - p->ioffset[1])/p->iscale[1];
	double nz = (z - p->ioffset[2])/p->iscale[2];
	double tmp[2];
	eval_nrpci(tmp, p, nx, ny, nz);
	result[0] = tmp[0] * p->scale[0] + p->offset[0];
	result[1] = tmp[1] * p->scale[1] + p->offset[1];
}



static double random_uniform(void)
{
	return rand()/(RAND_MAX+1.0);
}

static int main_trial(int c, char *v[])
{
	struct rpc p[1];
	nan_rpc(p);
	read_rpc_file_xml(p, "-");
	print_rpc(stderr, p, "p");
	for (int i = 0; i < 10; i++) {
		double x = random_uniform();
		double y = random_uniform();
		double z = random_uniform();
		double r[2], rr[2];
		eval_nrpc(r, p, x, y, z);
		eval_nrpci(rr, p, r[0], r[1], z);
		fprintf(stderr, "(%g %g %g) => (%g %g) => (%g %g)\n",
				x, y, z, r[0], r[1], rr[0], rr[1]);
		//fprintf(stderr, "%g\n", hypot(rr[0]-x, rr[1]-y));
	}
	return 0;
}

static int main_rpcline(int c, char *v[])
{
	if (c != 14 && c != 13) {
		fprintf(stderr, "usage:\n\t"
			"%s imga imgb a0x a0y b0x b0y rpca rpcb x y h0 hf [out]"
		//        0 1    2    3   4   5   6   7    8    9 10 11 12 13
			"\n", *v);
		return EXIT_FAILURE;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	int offset_a[2] = {atoi(v[3]), atoi(v[4])};
	int offset_b[2] = {atoi(v[5]), atoi(v[6])};
	char *filename_rpca = v[7];
	char *filename_rpcb = v[8];
	double basepoint[2] = {atof(v[9]), atof(v[10])};
	double hrange[2] = {atof(v[11]), atof(v[12])};

	struct rpc rpca[1]; read_rpc_file_xml(rpca, filename_rpca);
	struct rpc rpcb[1]; read_rpc_file_xml(rpcb, filename_rpcb);
	//print_rpc(stderr, rpca, "a");
	//print_rpc(stderr, rpcb, "b");
	int nh = 20;
	for (int i = 0; i < nh; i++)
	{
		double ix = basepoint[0];
		double iy = basepoint[1];
		double x = ix + offset_a[0];
		double y = iy + offset_a[1];
		double z = hrange[0] + i * (hrange[1] - hrange[0])/(nh-1);
		double r[2], rr[2];
		fprintf(stderr, "(%g %g %g) => ", ix, iy, z);
		eval_rpc(r, rpca, x, y, z);
		eval_rpci(rr, rpcb, r[0], r[1], z);
		double ox = rr[0] - offset_b[0];
		double oy = rr[1] - offset_b[1];
		//fprintf(stderr, "\n\tr[0] = %g\n", r[0]);
		//fprintf(stderr, "\tr[1] = %g\n", r[1]);
		//fprintf(stderr, "\trr[0] = %g\n", rr[0]);
		//fprintf(stderr, "\trr[1] = %g\n", rr[1]);
		fprintf(stderr, "(%g %g)\n", ox, oy);
	}
}

int main(int c, char *v[])
{
	return main_rpcline(c, v);
}
