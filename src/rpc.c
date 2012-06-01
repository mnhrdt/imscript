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
	FILE *f = xfopen(filename, "r");
	int n = 0x100, o = 1;
	while (1) {
		char line[n], tag[n], *sl = fgets(line, n, f);;
		if (!sl) break;
		if (o && 0 == strhas(line, "<Inverse_Model>")) o = 0;
		tag[0] = 'i';
		double x = get_xml_tagged_number(tag+1, line);
		if (isfinite(x)) {
			fprintf(stderr, "%s [%d]: %g\n", tag+o, o, x);
			add_tag_to_rpc(p, tag+o, x);
		}
	}
	xfclose(f);
}

#define FORI(n) for (int i = 0; i < (n); i++)

void print_rpc(FILE *f, struct rpc *p)
{
	FORI(20) fprintf(f, "numx[%d] = %g\n", i, p->numx[i]);
	FORI(20) fprintf(f, "denx[%d] = %g\n", i, p->denx[i]);
	FORI(20) fprintf(f, "numy[%d] = %g\n", i, p->numy[i]);
	FORI(20) fprintf(f, "deny[%d] = %g\n", i, p->deny[i]);
	FORI(20) fprintf(f, "inumx[%d] = %g\n", i, p->inumx[i]);
	FORI(20) fprintf(f, "idenx[%d] = %g\n", i, p->idenx[i]);
	FORI(20) fprintf(f, "inumy[%d] = %g\n", i, p->inumy[i]);
	FORI(20) fprintf(f, "ideny[%d] = %g\n", i, p->ideny[i]);
	FORI(3) fprintf(stderr, "scale[%d] = %g\n", i, p->scale[i]);
	FORI(3) fprintf(stderr, "offset[%d] = %g\n", i, p->offset[i]);
	FORI(3) fprintf(stderr, "iscale[%d] = %g\n", i, p->iscale[i]);
	FORI(3) fprintf(stderr, "ioffset[%d] = %g\n", i, p->ioffset[i]);

}

// evaluate a polynomial of degree 3
double eval_pol20(double c[20], double x, double y, double z)
{
	double r = 0;
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
}

// evaluate the normalized inverse rpc model
static void eval_nrpci(double *result,
		struct rpc *p, double x, double y, double z)
{
	double numx = eval_pol20(p->inumx, x, y, z);
	double denx = eval_pol20(p->idenx, x, y, z);
	double numy = eval_pol20(p->inumy, x, y, z);
	double deny = eval_pol20(p->ideny, x, y, z);
	result[0] = numx/denx;
	result[1] = numy/deny;
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
	result[1] = tmp[0] * p->iscale[1] + p->ioffset[1];
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
	result[1] = tmp[0] * p->scale[1] + p->offset[1];
}


// set all the values of an rpc model to NAN
static void nan_rpc(struct rpc *p)
{
	size_t nd = sizeof*p/sizeof(double);
	double *t = (double*)p;
	for (int i = 0; i < nd; i++)
		t[i] = NAN;
}

int main(int c, char *v[])
{
	struct rpc p[1];
	nan_rpc(p);
	read_rpc_file_xml(p, "-");
	print_rpc(stderr, p);
	return 0;
}
