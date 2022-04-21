// palette: a program to tranform from GRAY FLOATS to RGB BYTES
// SYNOPSIS:
// 	palette FROM TO PALSPEC [in.gray [out.rgb]]
//
//	FROM    = real value (-inf=min)
//	TO      = real value (inf=max)
//	PALSPEC = named palette (HOT, GLOBE, RELIEF, etc)j

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>


#define PALSAMPLES 1024
//#define PALSAMPLES 256

// Idea: The function
//
// 	palette : R -> rgb
//
// is described as a composition of two steps
//
// 	quantization : R -> {0, ..., PALSAMPLES-1}
// 	lookup : {0, ..., PALSAMPLES} -> rgb
//
// The first step is an affine function followed by a rounding, and the second
// step is the evaluation of a look-up table.  The proposed implementation is
// kept simple and the size of the look-up table is fixed.




#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"

struct palette {
	uint8_t t[3*PALSAMPLES];
	float m, M;
};

static void fprint_palette(FILE *f, struct palette *p)
{
	fprintf(f, "palette %g %g\n", p->m, p->M);
	for (int i = 0; i < PALSAMPLES; i++)
		fprintf(f, "\tpal[%d] = {%d %d %d}\n",
			i, p->t[3*i+0], p->t[3*i+1], p->t[3*i+2]);
}

static void fill_palette_gray(uint8_t *t)
{
	float factor = 255.0 / PALSAMPLES;
	for (int i = 0; i < PALSAMPLES; i++)
	for (int j = 0; j < 3; j++)
		t[3*i+j] = factor * i;
}

static void fill_palette_hot(uint8_t *t)
{
	float factor = 255.0 / PALSAMPLES;
	for (int i = 0; i < PALSAMPLES; i++)
	{
		t[3*i+0] = 255;
		if (i < 85)
			t[3*i+0] = 3 * factor * i;
		t[3*i+1] = factor * i;
		t[3*i+2] = factor * i;
	}
}

static void fill_palette_with_nodes(struct palette *p, float *n, int nn)
{
	if (nn < 2) fail("at least 2 nodes required");
	float factor = (PALSAMPLES-1.0)/(nn-1.0);
	//fprintf(stderr, "1pfnnm1 = %g\n", 1 + factor*(nn-1));
	assert(PALSAMPLES == round(1+factor*(nn-1)));
	for (int i = 0; i < nn-1; i++)
	{
		float posi = n[4*i];
		float posnext = n[4*(i+1)];
		//fprintf(stderr, "fwpn node %d : %g %g\n", i, posi, posnext);
		if (posi > posnext)
			fail("palette nodes must be in increasing order");
		int idxi = factor * i;
		int idxnext = factor * (i+1);
		assert(idxi <= idxnext);
		assert(idxnext < PALSAMPLES);
		for (int j = 0; j <= (idxnext-idxi); j++)
		{
			float a = j*1.0/(idxnext - idxi);
			p->t[3*(idxi+j)+0] = (1-a)*n[4*i+1] + a*n[4*(i+1)+1];
			p->t[3*(idxi+j)+1] = (1-a)*n[4*i+2] + a*n[4*(i+1)+2];
			p->t[3*(idxi+j)+2] = (1-a)*n[4*i+3] + a*n[4*(i+1)+3];
		}
	}
	p->m = n[0];
	p->M = n[4*(nn-1)];

	// dirty hack: copy the penultimate position of the palette to the last
	for (int j = 0; j < 3; j++)
		p->t[3*(PALSAMPLES-1)+j] = p->t[3*(PALSAMPLES-2)+j];
}

static void set_node_positions_linearly(float *n, int nn, float m, float M)
{
	float beta = (M - m)/(nn - 1);
	for (int i = 0; i < nn; i++)
		n[4*i] = m + beta * i;
}

static float nodes_cocoterrain[] = {
	1,   0,   0,   0, // black
	2, 255,   0, 255, // magenta
	3,   0,   0, 255, // blue
	4,   0, 255, 255, // cyan
	5,   0, 255,   0, // green
	6, 170,  84,   0, // brown
	7, 255, 255, 255, // white
};

static float nodes_dem[] = {
	0       ,  0       , 96.9994  , 70.9997  ,
	2.6010  , 16.0012  ,121.9997  , 46.9990 ,
	26.0100 , 231.9990 , 215.0007 , 125.0010,
	62.4495 , 160.9993 ,  67.0012 ,        0,
	88.4595 , 158.0006 ,        0 ,        0,
	145.7070,  109.9993,  109.9993,  109.9993,
	208.1565,  255.0000,  255.0000,  255.0000,
	255.0000,  255.0000,  255.0000,  255.0000,
};

static float nodes_nice[] = {
	1,   0,   0, 255, // blue
	2, 255, 255, 255, // white
	3, 255,   0,   0, // red
};

static float nodes_nnice[] = {
	1, 255,   0,   0, // red
	2,   0,   0,   0, // white
	3,   0, 255,   0, // blue
};

static float nodes_white[] = {
	1, 255, 255, 255,
	2, 255, 255, 255,
	3, 255, 255, 255,
};

static float nodes_black[] = {
	1, 0, 0, 0,
	2, 0, 0, 0,
	3, 0, 0, 0,
};

static float *get_gpl_nodes(char *filename, int *n)
{
	int bufsize = 0xff, nnodes = 0;
	char buf[bufsize];
	FILE *f = xfopen(filename, "r");
	static float nodes[4*PALSAMPLES];
	while(fgets(buf, bufsize, f) && nnodes < PALSAMPLES) {
		//fprintf(stderr, "s = \"%s\"\n", buf);
		float *t = nodes + 4*nnodes;
		*t = nnodes;
		if (3 == sscanf(buf, "%g %g %g", t+1, t+2, t+3)) {
			//fprintf(stderr, "\tnod[%d] = %g %g %g\n", nnodes, t[1], t[2], t[3]);
			nnodes += 1;
		}
	}
	xfclose(f);
	*n = nnodes;
	return nodes;
}

static float *get_gpf_nodes(char *filename, int *n)
{
	int bufsize = 0xff, nnodes = 0;
	char buf[bufsize];
	FILE *f = xfopen(filename, "r");
	static float nodes[4*PALSAMPLES];
	while(fgets(buf, bufsize, f) && nnodes < PALSAMPLES) {
		//fprintf(stderr, "s = \"%s\"\n", buf);
		float *t = nodes + 4*nnodes;
		*t = nnodes;
		if (4 == sscanf(buf, "%g %g %g %g", t, t+1, t+2, t+3)) {
			t[1] = round(t[1] * 255);
			t[2] = round(t[2] * 255);
			t[3] = round(t[3] * 255);
			//fprintf(stderr, "\tnod[%d]{%g} = %g %g %g\n", nnodes, t[0], t[1], t[2], t[3]);
			nnodes += 1;
		}
	}
	xfclose(f);
	*n = nnodes;
	return nodes;
}

static bool hassuffix(const char *s, const char *suf)
{
	int len_s = strlen(s);
	int len_suf = strlen(suf);
	if (len_s < len_suf)
		return false;
	return 0 == strcmp(suf, s + (len_s - len_suf));
}

static void fill_palette(struct palette *p, char *s, float m, float M)
{
	p->m = m;
	p->M = M;
	if (0 == strcmp(s, "gray"))
		fill_palette_gray(p->t);
	else if (0 == strcmp(s, "hot"))
		fill_palette_hot(p->t);
	else if (0 == strcmp(s, "cocoterrain")) {
		set_node_positions_linearly(nodes_cocoterrain, 7, m, M);
		fill_palette_with_nodes(p, nodes_cocoterrain, 7);
	} else if (0 == strcmp(s, "dem")) {
		set_node_positions_linearly(nodes_dem, 8, m, M);
		fill_palette_with_nodes(p, nodes_dem, 8);
	} else if (0 == strcmp(s, "nice")) {
		set_node_positions_linearly(nodes_nice, 3, m, M);
		fill_palette_with_nodes(p, nodes_nice, 3);
	} else if (0 == strcmp(s, "nnice")) {
		set_node_positions_linearly(nodes_nnice, 3, m, M);
		fill_palette_with_nodes(p, nodes_nnice, 3);
	} else if (0 == strcmp(s, "white")) {
		set_node_positions_linearly(nodes_white, 3, m, M);
		fill_palette_with_nodes(p, nodes_white, 3);
	} else if (0 == strcmp(s, "black")) {
		set_node_positions_linearly(nodes_black, 3, m, M);
		fill_palette_with_nodes(p, nodes_black, 3);
	} else if (hassuffix(s, ".gpl")) {
		int nnodes;
		float *nodes = get_gpl_nodes(s, &nnodes);
		set_node_positions_linearly(nodes, nnodes, m, M);
		fill_palette_with_nodes(p, nodes, nnodes);
	} else if (hassuffix(s, ".gpf")) {
		int nnodes;
		float *nodes = get_gpf_nodes(s, &nnodes);
		set_node_positions_linearly(nodes, nnodes, m, M);
		fill_palette_with_nodes(p, nodes, nnodes);
	}
	else fail("unrecognized palette \"%s\"", s);
}

static void get_palette_color(uint8_t *rgb, struct palette *p, float x)
{
	if (isnan(x)) {
		rgb[0] = rgb[1] = rgb[2] = 0;
		return;
	}
	if (!isfinite(x)) {
		rgb[0] = rgb[1] = rgb[2] = 255;
		return;
	}
	int ix = round((PALSAMPLES-1)*(x - p->m)/(p->M - p->m));
	if (ix < 0) ix = 0;
	if (ix >= PALSAMPLES) ix = PALSAMPLES - 1;
	rgb[0] = p->t[3*ix+0];
	rgb[1] = p->t[3*ix+1];
	rgb[2] = p->t[3*ix+2];
}

#include "smapa.h"
SMART_PARAMETER_SILENT(PALMAXEPS,0)

static void get_min_max(float *min, float *max, float *x, int n)
{
	float m = INFINITY, M = -m;
	for (int i = 0; i < n; i++)
		if (isfinite(x[i]))
		{
			m = fmin(m, x[i]);
			M = fmax(M, x[i]);
		}
	if (min) *min = m;
	if (max) *max = M+PALMAXEPS();
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static int get_two_quantiles(float *m, float *M, float *x, int n,
		float qm, float qM)
{
	float *tx = xmalloc(n*sizeof*tx);
	int N = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(x[i]))
			tx[N++] = x[i];
	if (N < 1) {
		fprintf(stderr, "too many NANs %d/%d", N, n);
		abort();
	}
	qsort(tx, N, sizeof*tx, compare_floats);
	int im = round(qm*N);
	int iM = round(qM*N);
	if (im < 0) im = 0;
	if (im >= N) im = N-1;
	if (iM < 0) iM = 0;
	if (iM >= N) iM = N-1;
	*m = tx[im];
	*M = tx[iM];
	free(tx);
	return N;
}
static int get_signed_quantile(float *t, float *x, int n, float p)
{
	float *tx = xmalloc(n*sizeof*tx);
	int N = 0;
	for (int i = 0; i < n; i++)
		if (!isnan(x[i]))
			tx[N++] = fabs(x[i]);
	if (N < 1) {
		fprintf(stderr, "too many NANs %d/%d", N, n);
		abort();
	}
	qsort(tx, N, sizeof*tx, compare_floats);
	int i = round(p*N);
	if (i < 0) i = 0;
	if (i >= N) i = N-1;
	*t = tx[i];
	free(tx);
	return N;
}

void parse_from_to(float *out_from, float *out_to, float *x, int n,
		char *from_id, char *to_id)
{
	char *mp, *Mp;
	float m = strtof(from_id, &mp);
	float M = strtof(to_id,   &Mp);
	fprintf(stderr, "parse m,M = %g %g\n", m, M);
	if (*mp != '%' && *Mp != '%')  // regular case, no percentiles
	{
		if (!isfinite(m)) get_min_max(&m, 0, x, n);
		if (!isfinite(M)) get_min_max(0, &M, x, n);
		*out_from = m;
		*out_to = M;
		return;
	}
	if (*mp == '%' && *Mp == '%') // two-sided percentile
	{
		get_two_quantiles(out_from, out_to, x, n, m/100, M/100);
		return;
	}
	if (m == 0 && *Mp == '%') // one-sided percentile: force signed palette
	{
		float t;
		get_signed_quantile(&t, x, n, M/100);
		*out_from = -t;
		*out_to = t;
		return;
	}
	if (M == 0 && *mp == '%') // one-sided percentile: force signed palette
	{
		float t;
		get_signed_quantile(&t, x, n, 1-m/100);
		*out_from = -t;
		*out_to = t;
		return;
	}
}

void apply_palette(uint8_t *y, float *x, int n, char *s, float m, float M)
{
	struct palette p[1];
	fill_palette(p, s, m, M);

	//fprint_palette(stderr, p);

	for (int i = 0; i < n; i++)
		get_palette_color(y + 3*i, p, x[i]);
}


#define PALETTE_MAIN

#ifdef PALETTE_MAIN
#include "iio.h"
#include "xmalloc.c"
#include "pickopt.c"

#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


SMART_PARAMETER_SILENT(PLEGEND_WIDTH,64)
SMART_PARAMETER_SILENT(PLEGEND_HEIGHT,256)
SMART_PARAMETER_SILENT(PLEGEND_MARGIN_LEFT,12)
SMART_PARAMETER_SILENT(PLEGEND_MARGIN_RIGHT,36)
SMART_PARAMETER_SILENT(PLEGEND_MARGIN_TOP,12)
SMART_PARAMETER_SILENT(PLEGEND_MARGIN_BOTTOM,12)
SMART_PARAMETER_SILENT(PLEGEND_TICKWIDTH,3)
SMART_PARAMETER_SILENT(PLEGEND_TEXT_XOFFSET,4)
SMART_PARAMETER_SILENT(PLEGEND_TEXT_YOFFSET,0)
SMART_PARAMETER_SILENT(PLEGEND_NTICKS,3)
SMART_PARAMETER_SILENT(PLEGEND_REVERSE,0)
SMART_PARAMETER_SILENT(PLEGEND_FACTOR,1)
SMART_PARAMETER_SILENT(PLEGEND_HIDE,0)
void save_legend(char *filename_legend, char *palette_id, float m, float M)
{
	m *= PLEGEND_FACTOR();
	M *= PLEGEND_FACTOR();

	// palette and font structs
	struct bitmap_font f[1] = {reformat_font(*xfont_7x14B, UNPACKED)};
	struct palette     p[1]; fill_palette(p, palette_id, m, M);

	// sizes, margins and positions
	int w = PLEGEND_WIDTH();
	int h = PLEGEND_HEIGHT();
	int m_l = PLEGEND_MARGIN_LEFT();
	int m_r = PLEGEND_MARGIN_RIGHT();
	int m_t = PLEGEND_MARGIN_TOP();
	int m_b = PLEGEND_MARGIN_BOTTOM();
	int p_i = m_l;     // left   boundary of colored part;
	int p_j = m_t;     // top    boundary of colored part
	int q_i = w - m_r; // right  boundary of colored part
	int q_j = h - m_b; // bottom boundary of colored part

	// image with the legend
	uint8_t *rgb = malloc(3*w*h);

	// fill background
	for (int i = 0; i < 3*w*h; i++)
		rgb[i] = 255;

	// transformation "x -> alpha * j + beta" from positions to values
	float alpha = (M - m) / (p_j - q_j);
	float beta  = m - alpha * q_j;
	if (PLEGEND_REVERSE()) {
		alpha = (m - M) / (p_j - q_j);
		beta = M - alpha * q_j;
	}

	// fill legend colors
	for (int j = p_j; j < q_j; j++)
	for (int i = p_i; i < q_i; i++)
	{
		//float x = m + ((M - m) * (j - p_j)) / (q_j - p_j);
		float x = alpha * j + beta;
		if (i == 64) fprintf(stderr, "j=%d x=%g\n", j, x);
		get_palette_color(rgb + 3*(j*w+i), p, x);
	}

	// border (1-pix black border)
	for (int l = 0; l < 3; l++) {
		for (int j = p_j; j < q_j; j++) {
			rgb[3*(j*w+p_i-1)+l] = 0;
			rgb[3*(j*w+q_i+0)+l] = 0;
		}
		for (int i = p_i-1; i < q_i; i++) {
			rgb[3*((p_j-1)*w+i)+l] = 0;
			rgb[3*((q_j+0)*w+i)+l] = 0;
		}
	}

	// ticks and numbers
	float ticks[][2] = {  // table with tick values and positions
		{M, p_j-1},
		{m, q_j},
		{(m+M)/2, (p_j+q_j)/2},
		{(3*m+M)/4, (3*p_j+q_j)/4},
		{(m+3*M)/4, (p_j+3*q_j)/4}
		// TODO: put more ticks (?)
	};
	int nticks = PLEGEND_NTICKS();
	if (nticks != 2 && nticks != 3 && nticks != 5)
		nticks = 2;
	for (int k = 0; k < nticks; k++)
	{
		float x = ticks[k][0]; // tick value
		int   j = ticks[k][1]; // tick position inside the legend
		if (PLEGEND_REVERSE()) x *= -1;

		// draw tick
		for (int i = q_i; i < q_i+1+PLEGEND_TICKWIDTH(); i++)
		for (int l = 0; l < 3; l++)
			rgb[3*(j*w+i)+l] = 0;

		// draw number associated to this tick
		char buf[0x100];
		uint8_t bg[3] = { 255, 255, 255}, fg[3] = {0, 0, 0};
		int X = PLEGEND_HIDE();
		if (!X)
			snprintf(buf, sizeof buf, "%g", x);
		else {
			snprintf(buf, sizeof buf, "XX%2.2lf", fmod(x,100));
		}
		int pos_i = q_i + PLEGEND_TICKWIDTH() + PLEGEND_TEXT_XOFFSET();
		int pos_j = j - f->height/2 + PLEGEND_TEXT_YOFFSET();
		put_string_in_rgb_image(rgb,w,h, pos_i, pos_j, fg,bg,0,f,buf);
	}

	// save legend into file
	iio_write_image_uint8_vec(filename_legend, rgb, w, h, 3);
	free(rgb);
}
//int main_palette(int c, char *v[])
//{
//	char *filename_legend = pick_option(&c, &v, "l", "");
//	if (c != 4 && c != 5 && c != 6 ) {
//		fprintf(stderr, "usage:\n\t%s from to pal [in [out]]\n", *v);
//		//                         0  1    2   3    4   5
//		return 1;
//	}
//	float from = atof(v[1]);
//	float to = atof(v[2]);
//	char *palette_id = v[3];
//	char *filename_in = c > 4 ? v[4] : "-";
//	char *filename_out = c > 5 ? v[5] : "-";
//
//	int w, h;
//	float *x = iio_read_image_float(filename_in, &w, &h);
//	uint8_t *y = xmalloc(3*w*h);
//
//	apply_palette(y, x, w*h, palette_id, &from, &to);
//
//	iio_write_image_uint8_vec(filename_out, y, w, h, 3);
//
//	if (*filename_legend)
//		save_legend(filename_legend, palette_id, from, to);
//
//	free(x);
//	free(y);
//	return 0;
//}
int main_palette2(int c, char *v[])
{
	char *filename_legend = pick_option(&c, &v, "l", "");
	if (c != 4 && c != 5 && c != 6 ) {
		fprintf(stderr, "usage:\n\t%s from to pal [in [out]]\n", *v);
		//                         0  1    2   3    4   5
		return 1;
	}
	char* from_id = v[1];
	char* to_id = v[2];
	char *palette_id = v[3];
	char *filename_in = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);
	uint8_t *y = xmalloc(3*w*h);

	float from, to;
	parse_from_to(&from, &to, x, w*h, from_id, to_id);
	fprintf(stderr, "from=%g to=%g\n", from, to);
	apply_palette(y, x, w*h, palette_id, from, to);

	iio_write_image_uint8_vec(filename_out, y, w, h, 3);

	if (*filename_legend)
		save_legend(filename_legend, palette_id, from, to);

	free(x);
	free(y);
	return 0;
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_palette2(c, v); }
#endif
#endif//PALETTE_MAIN
