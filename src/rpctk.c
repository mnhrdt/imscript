// toolkit for dealing with RPC functions
//
// rpctk localize rpc.xml < ijh.txt > llh.txt
// rpctk project  rpc.xml < llh.txt > ijh.txt
// rpctk fit              < xyhXY.txt > rpc.txt
// rpctk fitL             < ijhll.txt > rpc.xml
// rpctk fitP             < llhij.txt > rpc.xml
// rpctk fillL nx ny nh   > ijhll.txt
// rpctk fillP nx ny nh   > llhij.txt


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
	for (int j = 0; j < n; j++)
		fprintf(stderr, "ijhll[%d] = %g %g %g %g %g\n",
				j,
				ijhll[5*j+0],
				ijhll[5*j+1],
				ijhll[5*j+2],
				ijhll[5*j+3],
				ijhll[5*j+4]);
	double min_i   = min[0]; double max_i   = max[0];
	double min_j   = min[1]; double max_j   = max[1];
	double min_h   = min[2]; double max_h   = max[2];
	double min_lon = min[3]; double max_lon = max[3];
	double min_lat = min[4]; double max_lat = max[4];

fprintf(stderr, "mima_i = %g %g\n", min_i, max_i);
fprintf(stderr, "mima_j = %g %g\n", min_j, max_j);
fprintf(stderr, "mima_h = %g %g\n", min_h, max_h);
fprintf(stderr, "mima_lon = %g %g\n", min_lon, max_lon);
fprintf(stderr, "mima_lat = %g %g\n", min_lat, max_lat);

	// fill-in normalization factors
	r->scale[0]   = (max_i - min_i) / 2;
	r->scale[1]   = (max_j - min_j) / 2;
	r->scale[2]   = (max_h - min_h) / 2;
	r->offset[0]  = (max_i + min_i) / 2;
	r->offset[1]  = (max_j + min_j) / 2;
	r->offset[2]  = (max_h + min_h) / 2;
	r->iscale[0]  = (max_lon - min_lon) / 2;
	r->iscale[1]  = (max_lat - min_lat) / 2;
	r->iscale[2]  = (max_h   - min_h  ) / 2;
	r->ioffset[0] = (max_lon + min_lon) / 2;
	r->ioffset[1] = (max_lat + min_lat) / 2;
	r->ioffset[2] = (max_h   + min_h  ) / 2;

fprintf(stderr, "scale = %g %g %g\n", r->scale[0], r->scale[1], r->scale[2]);
fprintf(stderr, "offst = %g %g %g\n",r->offset[0],r->offset[1],r->offset[2]);
fprintf(stderr, "iscale = %g %g %g\n",r->iscale[0],r->iscale[1],r->iscale[2]);
fprintf(stderr,"ioffst = %g %g %g\n",r->ioffset[0],r->ioffset[1],r->ioffset[2]);

	// normalized input/outputs
	long double *ijh = xmalloc(3*n * sizeof*ijh);
	long double *lon = xmalloc(n * sizeof*lon);
	long double *lat = xmalloc(n * sizeof*lat);
	for (int i = 0; i < n; i++)
	{
		ijh[3*i+0] = (ijhll[5*i+0] - r->offset[0] ) / r->scale[0];
		ijh[3*i+1] = (ijhll[5*i+1] - r->offset[1] ) / r->scale[1];
		ijh[3*i+2] = (ijhll[5*i+2] - r->offset[2] ) / r->scale[2];
		lon[i]     = (ijhll[5*i+3] - r->ioffset[0]) / r->iscale[0];
		lat[i]     = (ijhll[5*i+4] - r->ioffset[1]) / r->iscale[1];
		if (n < 2000)
		{
			fprintf(stderr,"%d: %g %g %g %g %g -> "
					"%Lg %Lg %Lg %Lg %Lg\n",
					i,
					ijhll[5*i+0],
					ijhll[5*i+1],
					ijhll[5*i+2],
					ijhll[5*i+3],
					ijhll[5*i+4],
					ijh[3*i+0],
					ijh[3*i+1],
					ijh[3*i+2],
					lon[i],
					lat[i]);
		}
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
	fprintf(stderr, "e_lon = %g, e_lat = %g, e = %g\n", e_lon, e_lat, e);
	return e;
}

int main_rpctk_fillL(int c, char *v[])
{
	if (c != 5)
		return fprintf(stderr,"usage:\n\t"
				"%s rpc.txt nx ny nh > ijh.txt\n",*v);
	//                        0         1  2  3
	char *filename_rpc = v[1];
	int n[3] = { atoi(v[2]), atoi(v[3]), atoi(v[4]) };
	int nn = n[0] * n[1] * n[2];

	double *ijh = xmalloc(3 * nn * sizeof*ijh);

	struct rpc r[1];
	read_rpc_file_xml(r, filename_rpc);

	double Ai = (r->dmval[2] - r->dmval[0]) / (n[0] - 1);
	double Aj = (r->dmval[3] - r->dmval[1]) / (n[1] - 1);
	double Bi = r->dmval[0];
	double Bj = r->dmval[1];
	double Ak = 100;
	double Bk = 0;

	fprintf(stderr, "Ai Bi = %g %g\n", Ai, Bi);
	fprintf(stderr, "Aj Bj = %g %g\n", Aj, Bj);
	fprintf(stderr, "Ak Bk = %g %g\n", Ak, Bk);

	int cx = 0;
	for (int k = 0; k < n[2]; k++)
	for (int j = 0; j < n[1]; j++)
	for (int i = 0; i < n[0]; i++)
	{
		ijh[cx++] = Ai * i + Bi;
		ijh[cx++] = Aj * j + Bj;
		ijh[cx++] = Ak * k + Bk;
	}

	for (int i = 0; i < nn; i++)
		printf("%lf %lf %g\n", ijh[3*i+0], ijh[3*i+1], ijh[3*i+2]);

	free(ijh);
	return 0;
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

	free(ijhll);
	return 0;
}

#include <string.h>
int main(int c, char *v[])
{
	if (c < 2) goto end;
	if (0 == strcmp(v[1], "info"))     return main_rpctk_info(c-1, v+1);
	if (0 == strcmp(v[1], "localize")) return main_rpctk_localize(c-1, v+1);
	if (0 == strcmp(v[1], "project"))  return main_rpctk_project(c-1, v+1);
	if (0 == strcmp(v[1], "fitL"))     return main_rpctk_fitL(c-1, v+1);
	if (0 == strcmp(v[1], "fillL"))    return main_rpctk_fillL(c-1, v+1);
	//if (0 == strcmp(v[1], "fit"))      return main_rpctk_fit(c-1, v+1);
	//if (0 == strcmp(v[1], "fitP"))     return main_rpctk_fitP(c-1, v+1);
	//if (0 == strcmp(v[1], "fillP"))     return main_rpctk_filPL(c-1, v+1);
end:	return fprintf(stderr,
			"usage:\n\t%s {info|localize|project|fit} ...\n", *v);
}
