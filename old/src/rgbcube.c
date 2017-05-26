// simple color cube viewer
//
// usage:
// 	./rgbcube < image | display

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "iio.h"

#define XMALLOC_IS_MALLOC

#ifdef XMALLOC_IS_MALLOC
#define xmalloc malloc
#define xfree free
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


#ifndef FORI
#define FORI(n) for(int i=0;i<(int)(n);i++)
#endif//FORI

#ifndef FORJ
#define FORJ(n) for(int j=0;j<(int)(n);j++)
#endif//FORJ

#ifndef FORK
#define FORK(n) for(int k=0;k<(int)(n);k++)
#endif//FORK

#ifndef FORL
#define FORL(n) for(int l=0;l<(int)(n);l++)
#endif//FORL


// a smart parameter is just like a regular parameter, but it can be
// re-defined at the shell-environment.  Instead of
//
// 	#define NUMBER 42
// 	...
// 	printf("%g", NUMBER);
//
// do
// 	SMART_PARAMETER(NUMBER,42)
// 	...
// 	printf("%g", NUMBER());
//
// Notice that the environment only gets queried once, at the first use.
//
#define SMART_PARAMETER(n,v) static double n(void)\
{\
	static bool smapa_known_ ## n = false;\
	static double smapa_value_ ## n = v;\
	if (!smapa_known_ ## n)\
	{\
		fprintf(stderr,"scanning the environment for \"%s\"... ", #n);\
		int r;\
		char *sv = getenv(#n);\
		double y;\
		if (sv)\
			r = sscanf(sv, "%lf", &y);\
		if (sv && r == 1)\
		{\
			fprintf(stderr, "got value %g\n", y);\
			smapa_value_ ## n = y;\
		} else {\
			fprintf(stderr, "kept default value %g\n",\
					smapa_value_ ## n);\
		}\
		smapa_known_ ## n = true;\
	}\
	return smapa_value_ ## n;\
}

//static const char *emptystring = "";
//static const char *myname(void)
//{
//#define n 0x29a
//	static char buf[n];
//	pid_t p = getpid();
//	snprintf(buf, n, "/proc/%d/cmdline", p);
//	FILE *f = fopen(buf, "r");
//	if (!f) return emptystring;
//	int c, i = 0;
//	while ((c = fgetc(f)) != EOF && i < n) {
//#undef n
//		buf[i] = c ? c : ' ';
//		i += 1;
//	}
//	if (i) buf[i-1] = '\0';
//	fclose(f);
//	return buf;
//}

void error(const char *fmt, ...)
{
	va_list argp;
	//fprintf(stderr, "\nERROR(\"%s\"): ", myname());
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
//	if (global_hack_to_never_fail)
//	{
//		fprintf(stderr, "now wave a dead chicken and press enter\n");
//		getchar();
//		return;
//	}
//	exit(43);
#ifdef NDEBUG
	exit(-1);
#else
	//print_trace(stderr);
	exit(*(int *)0x43);
#endif
}

// assertion that all pointers are of the same size
// (so that some hacks below work)
static void assert_equal_pointers(void)
{
	size_t pvoid = sizeof(void *);
	size_t pchar = sizeof(char *);
	size_t pshort = sizeof(short *);
	size_t pint = sizeof(int *);
	size_t plong = sizeof(long *);
	size_t pllong = sizeof(long long *);
	size_t pfloat = sizeof(float *);
	size_t pdouble = sizeof(double *);
	size_t pldouble = sizeof(long double *);
	size_t pfunc = sizeof(&assert_equal_pointers);
	assert(pvoid == pchar);
	assert(pvoid == pshort);
	assert(pvoid == pint);
	assert(pvoid == plong);
	assert(pvoid == pllong);
	assert(pvoid == pfloat);
	assert(pvoid == pdouble);
	assert(pvoid == pldouble);
	assert(pvoid == pfunc);
}

// alloc a 2d matrix contiguously (wxh elements of size n)
void *matrix_build(int w, int h, size_t n)
{
	assert_equal_pointers();
	size_t p = sizeof(void *);
	char *r = xmalloc(h*p + w*h*n);
	for (int i = 0; i < h; i++)
		*(void **)(r + i*p) = r + h*p + i*w*n;
	return r;
}

// alloc a 3d matrix contiguously (wxhxd elements of size n)
void *matrix_build_3d(int w, int h, int d, size_t n)
{
	assert_equal_pointers();
	size_t p = sizeof(void *);
	void ***r = xmalloc(d*p + h*d*p + w*h*d*n);
	for (int i = 0; i < d; i++)
		r[i] = (void **)(r + d) + h*i;
	char *t = (char *)(r + d + h*d);
	for (int i = 0; i < d; i++)
		for (int j = 0; j < h; j++, t += w*n)
			r[i][j] = t;
	return r;
}



SMART_PARAMETER(CHISTO_BINS,64)
SMART_PARAMETER(CHISTO_NLOG,0)
SMART_PARAMETER(CHISTO_SMOOTH,0)
//SMART_PARAMETER(CHISTO_FULL,0)
//SMART_PARAMETER_SILENT(CHISTO_PROJ,0)
//SMART_PARAMETER(CHISTO_PROJ,0)
SMART_PARAMETER(CHISTO_RX,45)
SMART_PARAMETER(CHISTO_RY,60)
SMART_PARAMETER(CHISTO_RZ,72)
SMART_PARAMETER(CHISTO_VTRES,1)
SMART_PARAMETER(CHISTO_PD,6)
SMART_PARAMETER(CHISTO_PF,6)

static void matrix_product(float ab[3][3], float a[3][3], float b[3][3])
{
	FORI(3) FORJ(3) ab[i][j] = 0;
	FORI(3) FORJ(3) FORK(3) ab[i][j] += a[i][k] * b[k][j];
}
static void matrix_rotx(float a[3][3], float t)
{
	float ct = cos(t*M_PI/180.0);
	float st = sin(t*M_PI/180.0);
	FORI(3) FORJ(3) a[i][j] = 0;
	a[0][0] = 1;
	a[1][1] = a[2][2] = ct;
	a[1][2] = st;
	a[2][1] = -st;
}
static void matrix_roty(float a[3][3], float t)
{
	float ct = cos(t*M_PI/180.0);
	float st = sin(t*M_PI/180.0);
	FORI(3) FORJ(3) a[i][j] = 0;
	a[1][1] = 1;
	a[0][0] = a[2][2] = ct;
	a[0][2] = -st;
	a[2][0] = st;
}
static void matrix_rotz(float a[3][3], float t)
{
	float ct = cos(t*M_PI/180.0);
	float st = sin(t*M_PI/180.0);
	FORI(3) FORJ(3) a[i][j] = 0;
	a[2][2] = 1;
	a[1][1] = a[0][0] = ct;
	a[1][0] = st;
	a[0][1] = -st;
}
static void matrix_apply(float y[3], float a[3][3], float x[3])
{
	FORI(3) y[i] = 0;
	FORI(3) FORJ(3) y[i] += a[i][j] * x[j];
}

//static void mproj(int *o, float x[3], float r[3][3], int pside)
//{
//	float *flatr = (float *)r;
//	float y[3]; aplica_matriuf(y, flatr, x, 3);
//	FORI(2) y[i] = y[i] * CHISTO_PF()/(CHISTO_PD()+y[2]);
//	float py[2] = {y[0]/sqrt(3), y[1]/sqrt(3)};
//	//FORI(2) assert(py[i] >= -1 && py[i] <= 1);
//	FORI(2) py[i] = (py[i]+1)*0.5*pside;
//	FORI(2) o[i] = lrint(py[i]);
//	FORI(2) if (o[i] < 0) o[i] = 0;
//	FORI(2) if (o[i] >= pside) o[i] = pside-1;
//	if (0) {
//	fprintf(stderr, "mproj(%g %g %g)\n", x[0], x[1], x[2]);
//	fprintf(stderr, "\tby matrix (%g %g %g   %g %g %g    %g %g %g\n",
//			r[0][0], r[0][1], r[0][2],
//			r[1][0], r[1][1], r[1][2],
//			r[2][0], r[2][1], r[2][2]);
//	fprintf(stderr, "\tusing pside = %d\n", pside);
//	fprintf(stderr, "\ty = %g %g %g\n", y[0], y[1], y[2]);
//	fprintf(stderr, "\tpy = %g %g\n", py[0], py[1]);
//	fprintf(stderr, "\to = %d %d\n", o[0], o[1]);
//	}
//}

struct dot_to_display {
	float x, y, z;
	float r, g, b;
	float saturation;
};

static int compare_dots(const void *aa, const void *bb)
{
	const struct dot_to_display *a = (const struct dot_to_display *)aa;
	const struct dot_to_display *b = (const struct dot_to_display *)bb;
	int r = (a->z  > b->z) - (a->z < b->z);
	return -r;
}

struct list_of_dots {
	struct dot_to_display *t;
	int n, nmax;
};

static void start_list_of_dots(struct list_of_dots *l, int maxdot)
{
	l->t = xmalloc(maxdot * sizeof *l->t);
	l->n = 0;
	l->nmax = maxdot;
}

static void add_dot_to_list(struct list_of_dots *l,
		float x, float y, float z, float saturation)
{
	if (l->n >= l->nmax) error("out of dots");
	struct dot_to_display *d = l->t + l->n;
	d->r = x;
	d->g = y;
	d->b = z;
	d->saturation = saturation;
	l->n += 1;
}

// from rgb to xyz
static void project_dot(struct dot_to_display *d, float m[3][3],
		float view_f, float view_d, float pside)
{
	float x[3] = {d->r, d->g, d->b};
	FORI(3) x[i] = x[i] * 2.0/255.0 - 1;
	FORI(3) assert(x[i] >= -1);
	FORI(3) assert(x[i] <= 1);
	float y[3]; matrix_apply(y, m, x);
	FORI(2) y[i] = y[i] * view_f/(view_d+y[2]);
	float py[2] = {y[0]/sqrt(3), y[1]/sqrt(3)};
	FORI(2) py[i] = (py[i]+1)*0.5*pside;
	d->x = py[0];
	d->y = py[1];
	d->z = y[2];
}

//static void plot_dot(imatge_rgba *out, struct dot_to_display *d, int pside)
static void plot_dot(uint8_t (**out)[3], struct dot_to_display *d, int pside)
{
	//rgba_value mycol;
	//mycol.r = d->r;
	//mycol.g = d->g;
	//mycol.b = d->b;
	//mycol.a = 255;//d->b;
	int o[2] = {lrint(d->x), lrint(d->y)};
	FORI(2) if (o[i] < 0) o[i] = 0;
	FORI(2) if (o[i] >= pside) o[i] = pside-1;
	out[o[1]][o[0]][0] = d->r;
	out[o[1]][o[0]][1] = d->g;
	out[o[1]][o[0]][2] = d->b;
}

static void fill_histogram(float ***h, int bins, uint8_t (*x)[3], int nx)
{
	FORI(bins) FORJ(bins) FORK(bins) h[i][j][k] = 0;
	float binfac = 256.0/bins;
	for (int i = 0; i < nx; i++) {
		int r = x[i][0] / binfac;
		int g = x[i][1] / binfac;
		int b = x[i][2] / binfac;
		assert(r < bins && g < bins && b < bins);
		h[b][g][r] += 1;
	}

}

static void inplace_naive_smoothing(float ***x, int w, int h, int d)
{
	float ***th = matrix_build_3d(w, h, d, sizeof(float));
	FORK(d) FORJ(h) FORI(w)
	{
		int c = 0;
		float m = 0;
		for(int di = i-1; di <= i+1; di++)
		for(int dj = j-1; dj <= j+1; dj++)
		for(int dk = k-1; dk <= k+1; dk++)
		{
			if (di >= 0 && di < w)
			if (dj >= 0 && dj < h)
			if (dk >= 0 && dk < d)
			{
				m += x[dk][dj][di];
				c += 1;
			}
		}
		th[k][j][i] = m;
	}
	FORK(d) FORJ(h) FORI(w)
		x[k][j][i] = th[k][j][i];
	xfree(th);
}

static void beautify_histogram(float ***h, int bins)
{
	if (CHISTO_SMOOTH() > 0)
		inplace_naive_smoothing(h, bins, bins, bins);
	if (CHISTO_NLOG() > 0)
		FORI(bins) FORJ(bins) FORK(bins)
		       	h[i][j][k] = log(1 + h[i][j][k]);
}

// build a view of a 3d histogram
// input: histogram "h"
// output: image "out"
static void chistoview(uint8_t (**out)[3], int pside,
		float rx, float ry, float rz,
		float view_f, float view_d, float vtres,
		float ***h, int hside)
{
	//uint8_t (**out)[3] = matrix_build(pside, pside, 3);
	fprintf(stderr, "chistoview pside = %d, r = %g %g %g\n",pside,rx,ry,rz);
	fprintf(stderr, "chistoview vf vd vtres = %g %g %g\n",view_f,view_d,vtres);
	// number of points used to draw each edge of the cube
	int ndotsside = 300;

	// draw the background of the output image in gray
	FORJ(pside) FORI(pside) FORL(3)
		out[j][i][l] = 127;

	// count the total number of points
	int np = 0;
	FORK(hside) FORJ(hside) FORI(hside)
		if (h[k][j][i] >= vtres)
			np += 1;
	np += 12*ndotsside;

	struct list_of_dots l[1]; start_list_of_dots(l, np);
	int maxcolor = hside;
	int bw = 256 / maxcolor;

	// add the large enough bins
	FORK(hside) FORJ(hside) FORI(hside)
		if (h[k][j][i] >= vtres)
			add_dot_to_list(l, bw*i, bw*j, bw*k, h[k][j][i]);

	// add the edges of the cube
	FORI(ndotsside) {
		float p = pside - 1;
		float t = i * (p/(ndotsside-1.0));
		add_dot_to_list(l, 0, 0, t, INFINITY);
		add_dot_to_list(l, 0, t, 0, INFINITY);
		add_dot_to_list(l, t, 0, 0, INFINITY);
		add_dot_to_list(l, 0, p, t, INFINITY);
		add_dot_to_list(l, 0, t, p, INFINITY);
		add_dot_to_list(l, t, p, 0, INFINITY);
		add_dot_to_list(l, p, 0, t, INFINITY);
		add_dot_to_list(l, p, t, 0, INFINITY);
		add_dot_to_list(l, t, 0, p, INFINITY);
		add_dot_to_list(l, p, p, t, INFINITY);
		add_dot_to_list(l, p, t, p, INFINITY);
		add_dot_to_list(l, t, p, p, INFINITY);
	}

	// compute the approrpiate projection matrix
	float mrx[3][3]; matrix_rotx(mrx, rx);
	float mry[3][3]; matrix_roty(mry, ry);
	float mrz[3][3]; matrix_rotz(mrz, rz);
	float mrxy[3][3]; matrix_product(mrxy, mrx, mry);
	float m[3][3]; matrix_product(m, mrxy, mrz);

	// project all dots
	FORI(l->n)
		project_dot(l->t + i, m, view_f, view_d, pside);

	//sort by depth
	qsort(l->t, l->n, sizeof*l->t, compare_dots);

	// plot projected dots into output image
	FORI(l->n)
		plot_dot(out, l->t + i, pside);

	xfree(l->t);
}

static void draw_histogram(uint8_t (**d)[3], int dside, float ***h, int bins)
{
	chistoview(d, dside,
			CHISTO_RX(), CHISTO_RY(), CHISTO_RZ(),
			CHISTO_PF(), CHISTO_PD(), CHISTO_VTRES(),
			h, bins);
}

int main(void)
{
	int width, height, bins = CHISTO_BINS();
	uint8_t (*x)[3] = iio_read_image_uint8_rgb("-", &width, &height);
	fprintf(stderr, "got a %dx%d color image\n", width, height);

	float ***h = matrix_build_3d(bins, bins, bins, sizeof(float));
	fill_histogram(h, bins, x, width*height);
	beautify_histogram(h, bins);

	int dside = 256;
	uint8_t (**d)[3] = matrix_build(dside, dside, 3);
	draw_histogram(d, dside, h, bins);

	iio_write_image_uint8_matrix_rgb("-", d, dside, dside);

	xfree(d);
	xfree(h);
	xfree(x);
	return EXIT_SUCCESS;
}
