// toolkit for dealing with RPC functions
//
// rpctk localize rpc.xml < ijh.txt > llh.txt
// rpctk project  rpc.xml < llh.txt > ijh.txt
// rpctk fit              < xyhXY.txt > rpc.txt
// rpctk fitL             < ijhll.txt > rpc.xml
// rpctk fitP             < llhij.txt > rpc.xml


#include <stdio.h>

#include "xmalloc.c"

#include "rpcfit33.c"

// TODO: fix this shit
#define DONT_USE_TEST_MAIN
#include "ftr/rpc2.c"

int main_rpctk_info(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s file.rpc\n", *v);
	char *filename_rpc = v[1];
	struct rpc a[1];
	read_rpc_file_xml(a, filename_rpc);

	printf("RPC file \"%s\"\n", filename_rpc);
	printf("scale  = %g %g %g\n", a->scale[0], a->scale[1], a->scale[2]);
	printf("offset = %g %g %g\n", a->offset[0], a->offset[1], a->offset[2]);
	printf("numx   = %g %g ... %g\n", a->numx[0], a->numx[1], a->numx[19]);
	printf("denx   = %g %g ... %g\n", a->denx[0], a->denx[1], a->denx[19]);
	printf("numy   = %g %g ... %g\n", a->numy[0], a->numy[1], a->numy[19]);
	printf("deny   = %g %g ... %g\n", a->deny[0], a->deny[1], a->deny[19]);
	printf("iscale  = %g %g %g\n",a->iscale[0],a->iscale[1],a->iscale[2]);
	printf("ioffset = %g %g %g\n",*a->ioffset,a->ioffset[1],a->ioffset[2]);
	printf("inumx   = %g %g ... %g\n",a->inumx[0],a->inumx[1],a->inumx[19]);
	printf("idenx   = %g %g ... %g\n",a->idenx[0],a->idenx[1],a->idenx[19]);
	printf("inumy   = %g %g ... %g\n",a->inumy[0],a->inumy[1],a->inumy[19]);
	printf("ideny   = %g %g ... %g\n",a->ideny[0],a->ideny[1],a->ideny[19]);
	printf("dmval   = %g %g %g %g\n",
			a->dmval[0], a->dmval[1], a->dmval[2], a->dmval[3]);
	printf("imval  = %g %g %g %g\n",
			a->imval[0], a->imval[1], a->imval[2], a->imval[3]);
	return 0;
}

int main_rpctk_localize(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr,"usage:\n\t%s "
				"file.rpc < ijh.txt > llh.txt\n",*v);

	char *filename_rpc = v[1];
	struct rpc r[1];
	read_rpc_file_xml(r, filename_rpc);

	double ijh[3];
	while (3 == scanf("%lg %lg %lg\n", ijh, ijh+1, ijh+2))
	{
		double ll[2];
		rpc_localization(ll, r, ijh);
		printf("%lf %lf %lf\n", ll[0], ll[1], ijh[2]);
	}

	return 0;
}

int main_rpctk_project(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr,"usage:\n\t%s "
				"file.rpc < ijh.txt > llh.txt\n",*v);

	char *filename_rpc = v[1];
	struct rpc r[1];
	read_rpc_file_xml(r, filename_rpc);

	double ijh[3];
	while (3 == scanf("%lg %lg %lg\n", ijh, ijh+1, ijh+2))
	{
		double ll[2];
		rpc_projection(ll, r, ijh);
		printf("%lf %lf %lf\n", ll[0], ll[1], ijh[2]);
	}

	return 0;
}

// fit the localization function of an RPC given a (somewhat dense) list
// of localized pixels
static double rpc_fitL_full(struct rpc *r, double *ijhll, int n)
{
	// bounding boxes of input and output data
	double min[5], max[5];
	for (int i = 0; i < 5; i++)
	{
		min[i] =  INFINITY;
		max[i] = -INFINITY;
		for (int j = 0; j < n; j++)
			min[i] = fmin(min[i], ijhll[5*j+i]);
		for (int j = 0; j < n; j++)
			max[i] = fmax(max[i], ijhll[5*j+i]);
	}
	double min_i   = min[0], max_i   = max[0];
	double min_j   = min[1], max_j   = max[1];
	double min_h   = min[2], max_h   = max[2];
	double min_lon = min[3], max_lon = max[3];
	double min_lat = min[4], max_lat = max[4];

	// fill-in normalization factors
	r->scale[0] = 1;
	r->scale[1] = 1;
	r->scale[2] = 1;
	r->iscale[0] = 1;
	r->iscale[1] = 1;
	r->iscale[2] = 1;
	r->offset[0] = 1;
	r->offset[1] = 1;
	r->offset[2] = 1;
	r->ioffset[0] = 1;
	r->ioffset[1] = 1;
	r->ioffset[2] = 1;

	// normalized input/outputs
	long double *ijh = xmalloc(n * sizeof*ijh);
	long double *lon = xmalloc(n * sizeof*lon);
	long double *lat = xmalloc(n * sizeof*lat);
	for (int i = 0; i < n; i++)
	{
		ijh[3*i+0] = r->scale[0] * (ijhll[5*i+0] - r->offset[0]);
		ijh[3*i+1] = r->scale[1] * (ijhll[5*i+1] - r->offset[1]);
		ijh[3*i+2] = r->scale[2] * (ijhll[5*i+2] - r->offset[2]);
		lon[i]     = ijhll[5*i+3] / r->iscale[0] + r->ioffset[0];
		lat[i]     = ijhll[5*i+4] / r->iscale[1] + r->ioffset[1];
	}

	// fit the normalized model
	long double lon_p[20], lon_q[20], lat_p[20], lat_q[20];
	double e_lon = rpcfit33(lon_p, lon_q, ijh, lon, n);
	double e_lat = rpcfit33(lat_p, lat_q, ijh, lat, n);

	// cleanup
	free(ijh); free(lon); free(lat);

	// recover the coefficients of the normalized model
	for (int i = 0; i < 20; i++)
	{
		r->numx[i] = lon_p[i];
		r->denx[i] = lon_q[i];
		r->numy[i] = lat_p[i];
		r->deny[i] = lat_q[i];
	}

	// return the un-normalized error
	double e = e_lon * r->iscale[0] + e_lat * r->iscale[1];
	return e;
}

#include "parsenumbers.c"
int main_rpctk_fitL(int c, char *v[])
{
	if (c != 1)
		return fprintf(stderr,"usage:\n\t%s <ijhll.txt >rpc.txt\n", *v);

	int n;
	double *ijhll = read_ascii_doubles(stdin, &n);

	struct rpc r[1];
	nan_rpc(r);

	rpc_fitL_full(r, ijhll, n/5);

	print_rpc(stdout, r, "");

	return 0;
}

#include <string.h>
int main(int c, char *v[])
{
	if (c < 2) goto end;
	if (0 == strcmp(v[1], "info"))     return main_rpctk_info(c-1, v+1);
	if (0 == strcmp(v[1], "localize")) return main_rpctk_localize(c-1, v+1);
	if (0 == strcmp(v[1], "project"))  return main_rpctk_project(c-1, v+1);
	//if (0 == strcmp(v[1], "fit"))      return main_rpctk_fit(c-1, v+1);
	if (0 == strcmp(v[1], "fitL"))     return main_rpctk_fitL(c-1, v+1);
//	if (0 == strcmp(v[1], "fitP"))     return main_rpctk_fitP(c-1, v+1);
end:	return fprintf(stderr,
			"usage:\n\t%s {info|localize|project|fit} ...\n", *v);
}
