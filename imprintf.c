// usage:
// 	imprintf format [image]


#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"

#define xmalloc malloc

#define MAX_PIXELDIM IIO_MAX_DIMENSION

#define REQ_NOTHING 0
#define REQ_BASIC 1
#define REQ_SORTS 2
#define REQ_SORTP 4
#define REQ_SORTC 8
#define REQ_VNORM 16
#define REQ_SQUARES 32

struct conversion_specifier_and_its_data {
	char *name;
	char *shortname;

	int required_precomputation; // boolean mask of REQ_*

	int number_of_values;
	long double value[MAX_PIXELDIM];

	// only for conversion which take arguments:
	int argc;
	int argv[MAX_PIXELDIM];

	// flag for ignoring fields at various parts of the program
	bool selected;
};

struct printable_data {
	// the following table is initialized "statically" at the main function
	struct conversion_specifier_and_its_data *t;

	float *sorted_samples;
	float *sorted_vectors; // by comparing their norm
	float *sorted_colors;  // by comparing their samples
	float *squared_samples;
	float sum_squared_samples;

	//long double minsample, maxsample, avgsample, medsample;
	//long double minpixel[MAX_PIXELDIM];
	//long double maxpixel[MAX_PIXELDIM];
	//long double avgpixel[MAX_PIXELDIM];
	//long double medpixel[MAX_PIXELDIM];
	//long double nvectors, nscalars;

	char *string;
	//int strln;
	char *numberformat;
	char *vectorspacing;

	int action_table[0x100];
	int compuflag;

	float *x;
	int w, h, pd;
};

static int getidx(struct printable_data *p, char *id)
{
	int i = 0;
	while(p->t[i].name)
	{
		if (0 == strcmp(id, p->t[i].name))
			return i;
		i = i + 1;
	}
	exit(fprintf(stderr, "ERROR: bad %%id = \"%s\"\n", id));
}

static void setnumber(struct printable_data *p, char *id, long double x)
{
	int idx = getidx(p, id);
	p->t[idx].number_of_values = 1;
	p->t[idx].value[0] = x;
	//fprintf(stderr, "SETNUMBER[%d \"%s\"] = %g\n", idx, id, (double)x);
}

static void setnumbers(struct printable_data *p, char *id,long double *x, int n)
{
	int idx = getidx(p, id);
	assert(0 <= n && n <= MAX_PIXELDIM);
	p->t[idx].number_of_values = n;
	for (int i = 0; i < n; i++)
		p->t[idx].value[i] = x[i];
}

static int getnumbers(long double *x, struct printable_data *p, char *id)
{
	int idx = getidx(p, id);
	int r = p->t[idx].number_of_values;
	assert(r > 0);
	for (int i = 0; i < r; i++)
		x[i] = p->t[idx].value[i];
	return r;
}

static long double getnumber(struct printable_data *p, char *id)
{
	int idx = getidx(p, id);
	int n = p->t[idx].number_of_values;
	assert(n == 1);
	return p->t[idx].value[0];
}

static void config_printable_data(struct printable_data *p, char *fmt)
{
	// hack
	p->string = fmt;

	// more hack
	int idx = 0;
	for (int i = 0; i < 0x100; i++)
		p->action_table[i] = -1;
	while (p->t[idx].name) {
		int act = ((unsigned int)p->t[idx].shortname[0]) % 0x100;
		p->action_table[act] = idx;
		idx += 1;
	}

	// parse
	while(*fmt) {
		int c = *fmt++;
		if (c == '%') {
			c = *fmt++;
			idx = p->action_table[c%0x100];
			if (idx >= 0) {
				//printf("GOT: \"%s\"\n", p->t[idx].name);
				p->t[idx].selected = true;
			}
		}// else
		//	printf("PUT: '%c'\n", c);
	}
}

static float vnormf(float *x, int n)
{
	float r = 0;
	for (int i = 0; i < n; i++)
		r = hypot(r, x[i]);
	return r;
}

static void compute_stuff_nothing(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	setnumber(p, "width", w);
	setnumber(p, "height", h);
	setnumber(p, "depth", 1);
	setnumber(p, "pixeldim", pd);
	setnumber(p, "nsamples", w*h*pd);
	setnumber(p, "npixels", w*h);
	setnumber(p, "dimension", 2);
}

static void compute_stuff_basic(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	// sample basic stuff
	int ns = w * h * pd;
	float minsample = INFINITY, maxsample = -INFINITY;
	long double avgsample = 0;
	for (int i = 0; i < ns; i++)
	{
		float y = x[i];
		if (y < minsample) minsample = y;
		if (y > maxsample) maxsample = y;
		avgsample += y;
	}
	avgsample /= ns;
	setnumber(p, "minsample", minsample);
	setnumber(p, "maxsample", maxsample);
	setnumber(p, "avfsample", avgsample);


	// pixel basic stuff
	int np = w * h, minpixel_idx, maxpixel_idx;
	float minpixel = INFINITY, maxpixel = -INFINITY;
	long double avgpixel[MAX_PIXELDIM] = {0};
	int minidx, maxidx;
	for (int j = 0; j < pd; j++)
		avgpixel[j] = 0;
	for (int i = 0; i < np; i++)
	{
		float xnorm = vnormf(x + pd*i, pd);
		if (xnorm < minpixel) { minidx = i; minpixel = xnorm; }
		if (xnorm > maxpixel) { maxidx = i; maxpixel = xnorm; }
		for (int j = 0; j < pd; j++)
			avgpixel[j] += x[pd*i+j];
	}
	long double mipi[pd], mapi[pd];
	for (int j = 0; j < pd; j++) {
		mipi[j] = x[minidx*pd+j];
		mapi[j] = x[minidx*pd+j];
		avgpixel[j] /= np;
	}
	setnumbers(p, "minpixel", mipi, pd);
	setnumbers(p, "maxpixel", mapi, pd);
	setnumbers(p, "avgpixel", avgpixel, pd);
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void compute_stuff_sorts(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	int ns = w * h * pd;
	p->sorted_samples = xmalloc(ns * sizeof*p->sorted_samples);
	for (int i = 0; i < ns; i++)
		p->sorted_samples[i] = x[i];
	qsort(p->sorted_samples, ns, sizeof*p->sorted_samples, compare_floats);
}

static void compute_stuff_sortp(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	exit(fprintf(stderr, "ERROR: sortp not yet implemented\n"));
}

static void compute_stuff_sortc(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	exit(fprintf(stderr, "ERROR: sortc not yet implemented\n"));
}

static void compute_stuff_vnorm(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	exit(fprintf(stderr, "ERROR: vnorm not yet implemented\n"));
}

static void compute_stuff_squares(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	int n = w * h * pd;
	p->squared_samples = xmalloc(n * sizeof*p->squared_samples);
	long double hyp = 0;
	for (int i = 0; i < n; i++)
	{
		float sq = x[i] * x[i];
		p->squared_samples[i] = sq;
		hyp = hypotl(hyp, sq);
	}
	p->sum_squared_samples = hyp * hyp;
}


#define GETBIT(x,i) ((x)&(1<<(i)))

static void compute_printable_data(struct printable_data *p,
		float *x, int w, int h, int pd)
{
	p->x = x; p->w = w; p->h = h; p->pd = pd;
	int action_table[0x100], idx = 0;
	p->compuflag;
	for (int i = 0; i < 0x100; i++)
		action_table[i] = -1;
	while (p->t[idx].name) {
		if (p->t[idx].selected) {
			p->compuflag |= p->t[idx].required_precomputation;
			//fprintf(stderr, "going to compute \"%s\"\n",
			//		p->t[idx].name);
		}
		idx += 1;
	}

	//fprintf(stderr, "            ba ss sp sc vn sq\n");
	//fprintf(stderr, "compuflag = %d  %d  %d  %d  %d  %d\n",
	//		(bool)GETBIT(p->compuflag,0),
	//		(bool)GETBIT(p->compuflag,1),
	//		(bool)GETBIT(p->compuflag,2),
	//		(bool)GETBIT(p->compuflag,3),
	//		(bool)GETBIT(p->compuflag,4),
	//		(bool)GETBIT(p->compuflag,5)
	//       );

	compute_stuff_nothing(p, x, w, h, pd);
	if (REQ_BASIC & p->compuflag) compute_stuff_basic(p, x, w, h, pd);
	if (REQ_SORTS & p->compuflag) compute_stuff_sorts(p, x, w, h, pd);
	if (REQ_SORTP & p->compuflag) compute_stuff_sortp(p, x, w, h, pd);
	if (REQ_SORTC & p->compuflag) compute_stuff_sortc(p, x, w, h, pd);
	if (REQ_VNORM & p->compuflag) compute_stuff_vnorm(p, x, w, h, pd);
	if (REQ_SQUARES & p->compuflag) compute_stuff_squares(p, x, w, h, pd);
}

static void print_printable_datum(FILE *f, struct printable_data *p, int idx)
{
	//fprintf(f, "PRINTABLE DATUM[%d \"%s\"] = \"", idx, p->t[idx].name);
	int n = p->t[idx].number_of_values;
	for (int i = 0; i < n; i++) {
		fprintf(f, p->numberformat, (double)(p->t[idx].value[i]));
		if (i > 0)
			fprintf(f, "%s", p->vectorspacing);
	}
	//fprintf(f, "\"\n");
}

static void my_putchar(FILE *f, int c)
{
	fputc(c, f);
}

static void print_printable_data(FILE *f, struct printable_data *p)
{
	char *fmt = p->string;

	// parse
	while(*fmt) {
		int c = *fmt++;
		if (c == '%') {
			c = *fmt++; if (!c) break;
			int idx = p->action_table[c%0x100];
			if (idx >= 0) {
				//printf("_GOT: \"%s\"\n", p->t[idx].name);
				print_printable_datum(f, p, idx);
				p->t[idx].selected = true;
			}
		} else if (c == '\\') {
			c = *fmt++; if (!c) break;
			if (c == 'n') my_putchar(f, '\n');
			if (c == 't') my_putchar(f, '\t');
			if (c == '\\') my_putchar(f, '\\');
		} else
			my_putchar(f, c);
	}
}

static void imprintf_2d(FILE *f, char *fmt, float *x, int w, int h, int pd)
{
	struct printable_data p[1] = {{
		.t = (struct conversion_specifier_and_its_data []){
			{"width",      "w", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"height",     "h", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"depth",      "d", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"pixeldim",   "c", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"nsamples",   "n", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"npixels",    "N", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"getsample",  "p", REQ_NOTHING, 1, {0},-1, {0}, false},
			{"getpixel",   "P", REQ_NOTHING,-1, {0},-1, {0}, false},
			{"dimension",  "D", REQ_NOTHING, 1, {0}, 0, {0}, false},
			{"size",       "z", REQ_NOTHING, 1, {0}, 1, {0}, false},
			{"sizes",      "Z", REQ_NOTHING,-1, {0}, 1, {0}, false},
			{"minsample",  "i", REQ_BASIC,   1, {0}, 0, {0}, false},
			{"maxsample",  "a", REQ_BASIC,   1, {0}, 0, {0}, false},
			{"avgsample",  "v", REQ_BASIC,   1, {0}, 0, {0}, false},
			{"medsample",  "m", REQ_SORTS,   1, {0}, 0, {0}, false},
			{"percentile", "q", REQ_SORTS,   1, {0}, 1, {0}, false},
			{"minpixel",   "I", REQ_BASIC,  -1, {0}, 0, {0}, false},
			{"avgpixel",   "V", REQ_BASIC,  -1, {0}, 0, {0}, false},
			{"medpixel",   "M", REQ_SORTP,  -1, {0}, 0, {0}, false},
			{"Percentile", "Q", REQ_SORTP,   1, {0}, 1, {0}, false},
			{"nscalars",   "k", REQ_SORTS,   1, {0}, 0, {0}, false},
			{"nvectors",   "K", REQ_SORTC,   1, {0}, 0, {0}, false},
			{"rms",        "r", REQ_SQUARES, 1, {0}, 0, {0}, false},
			{NULL, NULL, 0, 0, {0}, 0, {0}, false} },
		.sorted_samples = NULL,
		.sorted_vectors = NULL,
		.sorted_colors = NULL,
		.squared_samples = NULL,
		//.nvectors = NAN,
		//.nscalars = NAN,
		.string = NULL,
		.numberformat = "%g",
		.vectorspacing = ", ",
		//.strln = 0,
	}};

	config_printable_data(p, fmt);
	compute_printable_data(p, x, w, h, pd);
	print_printable_data(f, p);
}

int main(int c, char *v[])
{
	if (c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s format [image]\n", *v);
		//                          0 1       2
		return EXIT_FAILURE;
	}
	char *format = v[1];
	char *finame = c > 2 ? v[2] : "-";
	int w, h, pd;
	float *x = iio_read_image_float_vec(finame, &w, &h, &pd);
	imprintf_2d(stdout, format, x, w, h, pd);
	return EXIT_SUCCESS;
}
