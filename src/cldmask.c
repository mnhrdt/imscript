//////////////////
/// INTERFACE ////
//////////////////


struct cloud_polygon {
	int n;                       // number of vertices
	double *v;                   // vertex coordinates (array of length 2*n)
};

struct cloud_mask {
	int n;                       // number of polygons
	struct cloud_polygon *t;     // array of polygons
	double low[2], up[2]; // "rectangle"
};

int read_cloud_mask_from_gml_file(struct cloud_mask *m, char *filename);

void clouds_mask_fill(float *img, int w, int h, struct cloud_mask *m);





///////////////////////
/// IMPLEMENTATION ////
///////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"

// read stream until character "stop" is found
// if EOF is reached, return NULL
// otherwise return line
static char *fgets_until(char *line, int n, FILE *f, int stop)
{
	int i = 0;
	while(1) {
		if (i >= n-1) break;
		int c = fgetc(f);
		if (c == EOF) return NULL;
		line[i] = c;
		if (c == stop) break;
		i += 1;
	}
	line[i+1] = '\0';
	//fprintf(stderr, "FGETS UNTIL %d \"%s\"\n", i, line);
	return line;
}

static void read_until_newline(FILE *f)
{
	while (1) {
		int c = fgetc(f);
		if (c == EOF || c == '\n')
			return;
	}
}

// like strcmp, but finds a needle
static int strhas(char *haystack, char *needle)
{
	char *r = strstr(haystack, needle);
	return r ? 0 : 1;
}

static void cloud_add_polygon(struct cloud_mask *m, struct cloud_polygon p)
{
	m->n += 1;
	m->t = xrealloc(m->t, m->n * sizeof*m->t);
	m->t[m->n-1] = p;
}

// needed only for checking consistency at the end
void free_cloud(struct cloud_mask *m)
{
	for (int i = 0; i < m->n; i++)
		free(m->t[i].v);
	free(m->t);
}

int read_cloud_mask_from_gml_file(struct cloud_mask *m, char *filename)
{
	m->n = 0;
	m->t = NULL;

	FILE *f = xfopen(filename, "r");
	int n = 0x100, nf;
	while(1) {
		char line[n], *sl = fgets_until(line, n, f, '>');
		if (!sl) break;
		if (0 == strhas(line, "lowerCorner")) {
			double *ff = read_ascii_doubles(f, &nf);
			if (nf == 2) for (int i = 0; i < 2; i++)
				m->low[i] = ff[i];
			free(ff);
		}
		if (0 == strhas(line, "upperCorner")) {
			double *ff = read_ascii_doubles(f, &nf);
			if (nf == 2) for (int i = 0; i < 2; i++)
				m->up[i] = ff[i];
			free(ff);
		}
		if (0 == strhas(line, "posList")) {
			struct cloud_polygon p;
			p.v = read_ascii_doubles(f, &p.n);
			p.n /= 2;
			cloud_add_polygon(m, p);
		}
		read_until_newline(f);
	}
	return 0;
}

static void putpixel_1(float *img, int w, int h, float x, float y, float v)
{
	int i = round(x);
	int j = round(y);
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	img[j*w+i] = v;
}


struct plot_data { float *x; int w, h; float c; };
static void plot_pixel(int i, int j, void *e)
{
	struct plot_data *d = e;
	putpixel_1(d->x, d->w, d->h, i, j, d->c);
}

#include "drawsegment.c"
static void plot_segment_gray(float *img, int w, int h,
		float a[2], float b[2], float c)
{
	int p[2] = {round(a[0]), round(a[1])};
	int q[2] = {round(b[0]), round(b[1])};
	struct plot_data d = { img, w, h, c };
	traverse_segment(p[0], p[1], q[0], q[1], plot_pixel, &d);
}

void clouds_mask_fill(float *img, int w, int h, struct cloud_mask *m)
{
	//fprintf(stderr, "CLOUD RECTANGLE %g %g %g %g\n", m->low[0], m->low[1], m->up[0], m->up[1]);
	for (int i = 0; i < m->n; i++)
	{
		struct cloud_polygon *p = m->t + i;
		//fprintf(stderr, "CLOUD %d/%d : %d :", i+1, m->n+1, p->n);
		for (int j = 0; j < p->n - 1; j++)
		{
			float a[2] = {p->v[2*j+0], p->v[2*j+1]};
			float b[2] = {p->v[2*j+2], p->v[2*j+3]};
			float A[2] = {
				(a[0] - m->low[0]) / (m->up[0] - m->low[0]) * w,
				(a[1] - m->low[1]) / (m->up[1] - m->low[1]) * h
			};
			float B[2] = {
				(b[0] - m->low[0]) / (m->up[0] - m->low[0]) * w,
				(b[1] - m->low[1]) / (m->up[1] - m->low[1]) * h
			};
			plot_segment_gray(img, w, h, A, B, 128);
			putpixel_1(img, w, h, A[0], A[1], 255);
			putpixel_1(img, w, h, B[0], B[1], 255);
			//fprintf(stderr, " (%g %g)", p->v[2*j], p->v[2*j+1]);
		}
		//fprintf(stderr, "\n");
	}
}


#include "iio.h"
int main(int c, char *v[])
{
	// read input arguments
	if (c != 4) {
		return fprintf(stderr, "usage:\n\t"
				"%s CLD_.clg factor image.png\n", *v);
		//                0 1        2      3
	}
	char *filename_clg = v[1];
	float factor = atof(v[2]);
	char *filename_out = v[3];

	// read input cloud file
	struct cloud_mask m[1];
	read_cloud_mask_from_gml_file(m, filename_clg);

	// acquire space for output image
	int w = 1000;
	int h = 1000;
	float *x = xmalloc(w*h*sizeof*x);
	for (int i = 0; i < w*h; i++)
		x[i] = 0;

	// draw masks over output image
	clouds_mask_fill(x, w, h, m);

	// save output image
	iio_save_image_float(filename_out, x, w, h);

	
	//cleanup
	free(x);
	free_cloud(m);
	return 0;
}
