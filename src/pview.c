// display points, pairs or triplets
// (optionally transformed, optionally with color masks)
// all programs below produce a color image as output
// the last two arguments "w h" are the size of the desired image
//
// pview points w h < points.txt
// pview pairs h1 ... h9 w h [mask.txt] < pairs.txt
// pview triplets w h < triplets.txt
// pview triplets w h [mask.txt] < triplets.txt
// pview epipolar f1 ... f9 w h < pairs.txt
// pview epipolar f1 ... f9 w h [mask.txt] < pairs.txt
//
// points.txt   = file with two columns of numbers (list of 2D points)
// pairs.txt    = file with four columns of numbers (list of 2D point pairs)
// triplets.txt = file with six columns of numbers (list of 2D point triplets)


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"
#include "drawsegment.c"

struct rgba_value {
	uint8_t r, g, b, a;
};
#define RGBA_BLACK   (struct rgba_value){.r=0x00,.g=0x00,.b=0x00,.a=0xff}
#define RGBA_RED     (struct rgba_value){.r=0xff,.g=0x00,.b=0x00,.a=0xff}
#define RGBA_GREEN   (struct rgba_value){.r=0x00,.g=0xff,.b=0x00,.a=0xff}
#define RGBA_BLUE    (struct rgba_value){.r=0x00,.g=0x00,.b=0xff,.a=0xff}
#define RGBA_MAGENTA (struct rgba_value){.r=0xff,.g=0x00,.b=0xff,.a=0xff}
#define RGBA_YELLOW  (struct rgba_value){.r=0xff,.g=0xff,.b=0x00,.a=0xff}
#define RGBA_PHANTOM (struct rgba_value){.r=30,.g=20,.b=10,.a=0xff}
#define RGBA_BRIGHT  (struct rgba_value){.r=129,.g=86,.b=43,.a=0xff}
#define RGBA_DARKGRAY (struct rgba_value){.r=0x30,.g=0x30,.b=0x30,.a=0xff}
#define RGBA_GRAY10 (struct rgba_value){.r=0x10,.g=0x10,.b=0x10,.a=0xff}
#define RGBA_GRAY20 (struct rgba_value){.r=0x20,.g=0x20,.b=0x20,.a=0xff}
#define RGBA_GRAY30 (struct rgba_value){.r=0x30,.g=0x30,.b=0x30,.a=0xff}
#define RGBA_GRAY40 (struct rgba_value){.r=0x40,.g=0x40,.b=0x40,.a=0xff}
#define RGBA_GRAY50 (struct rgba_value){.r=0x50,.g=0x50,.b=0x50,.a=0xff}



#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif//FORI

static bool inner_point(int w, int h, int x, int y)
{
	return (x>=0) && (y>=0) && (x<w) && (y<h);
}

static int main_viewp(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s sx sy < pairs.txt\n", *v);
		//                         0  1  2
		return EXIT_FAILURE;
	}
	int n;
	float *t = read_ascii_floats(stdin, &n);
	n /= 2;
	int sizex = atoi(v[1]);
	int sizey = atoi(v[2]);
	struct rgba_value (*x)[sizex] = xmalloc(4*sizex*sizey);
	FORI(sizex*sizey)
		x[0][i] = RGBA_BLACK;
	FORI(n)
	{
		int a = t[2*i+0];
		int b = t[2*i+1];
		if (inner_point(sizex, sizey, a, b))
			x[b][a] = RGBA_GREEN;
	}
	iio_save_image_uint8_vec("-", (uint8_t*)x, sizex, sizey, 4);
	free(t); free(x);
	return EXIT_SUCCESS;
}

void projective_map(double y[2], double A[9], double x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
	double y2 = A[6]*x[0] + A[7]*x[1] + A[8];
	y[0] /= y2;
	y[1] /= y2;
}

static void draw_phantom_pixel(int a, int b, void *pp)
{
	struct {int w, h; struct rgba_value *x;} *p = pp;
	if (inner_point(p->w, p->h, a, b))
		p->x[b*p->w+a] = RGBA_PHANTOM;
}

static void draw_brighter_pixel(int a, int b, void *pp)
{
	struct {int w, h; struct rgba_value *x;} *p = pp;
	if (inner_point(p->w, p->h, a, b))
		p->x[b*p->w+a] = RGBA_BRIGHT;
}


int main_viewpairs(int c, char *v[])
{
	if (c != 12 && c != 13) {
		fprintf(stderr, "usage:\n\t"
			"%s a b r c d s p q 1 sx sy [mask] < pairs.txt\n", *v);
		//       0  1 2 3 4 5 6 7 8 9 10 11 12
		return EXIT_FAILURE;
	}
	double A[9]; FORI(9) A[i] = atof(v[1+i]);
	int sizex = atoi(v[10]);
	int sizey = atoi(v[11]);
	int n;
	double (*p)[4] = (void*)read_ascii_doubles(stdin, &n);
	n /= 4;
	struct rgba_value (*o)[sizex] = xmalloc(sizex*sizey*4);
	struct {int w, h; struct rgba_value *x; } usr = {sizex, sizey, o[0]};
	bool bmask[n]; FORI(n) bmask[i] = c<13;
	if (c>12) { FILE *f = xfopen(v[12], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=t;}
		xfclose(f);
	}
	FORI(sizex*sizey) o[0][i] = RGBA_BLACK;
	FORI(n) {
		double *pxi = p[i];
		double *pyi = p[i]+2;
		double tt[2]; projective_map(tt, A, pxi);
		int t[2] = {tt[0], tt[1]};
		int z[2] = {pyi[0], pyi[1]};
		if (inner_point(sizex, sizey, t[0], t[1]) &&
				inner_point(sizex, sizey, z[0], z[1]))
			if (!bmask[i]) {
				traverse_segment(t[0], t[1], z[0], z[1],
						draw_phantom_pixel, &usr);
				o[t[1]][t[0]] = RGBA_RED;
				o[z[1]][z[0]] = RGBA_BLUE;
			}
	}
	FORI(n) {
		double tt[2]; projective_map(tt, A, p[i]);
		int t[2] = {tt[0], tt[1]};
		int z[2] = {p[i][2], p[i][3]};
		if (inner_point(sizex, sizey, t[0], t[1]) &&
				inner_point(sizex, sizey, z[0], z[1])
		   ) {
			if (bmask[i]) {
				traverse_segment(t[0], t[1], z[0], z[1],
						draw_brighter_pixel, &usr);
				o[t[1]][t[0]] = RGBA_RED;
				o[z[1]][z[0]] = RGBA_GREEN;
			}
		}
	}
	iio_save_image_uint8_vec("-", (uint8_t*)o, sizex, sizey, 4);
	return EXIT_SUCCESS;
}

int main_viewtrips(int c, char *v[])
{
	if (c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s sx sy [mask] < pairs.txt\n", *v);
		//                         0  1  2   3
		return EXIT_FAILURE;
	}
	int s[2], n; FORI(2) s[i] = atoi(v[1+i]);
	double (*p)[6] = (void*)read_ascii_doubles(stdin, &n);
	n /= 6;
	struct rgba_value (*o)[s[0]] = xmalloc(s[0]*s[1]*4);
	struct {int w, h; struct rgba_value *x; } usr = {s[0], s[1], o[0]};
	bool bmask[n]; FORI(n) bmask[i] = c<4;
	if (c>3) { FILE *f = xfopen(v[3], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=t;}
		xfclose(f);
	}
	FORI(s[0]*s[1]) o[0][i] = RGBA_BLACK;
	FORI(n) {
		int t[2] = {p[i][0], p[i][1]};
		int z[2] = {p[i][2], p[i][3]};
		int Z[2] = {p[i][4], p[i][5]};
		if (inner_point(s[0], s[1], t[0], t[1]) &&
				inner_point(s[0], s[1], z[0], z[1]) &&
				inner_point(s[0], s[1], Z[0], Z[1]))
		{
			void (*pix)(int,int,void*) = (bmask[i]&&c>3)?
				draw_brighter_pixel:draw_phantom_pixel;
			traverse_segment(t[0], t[1], z[0], z[1],
					pix, &usr);
			traverse_segment(z[0], z[1], Z[0], Z[1],
					pix, &usr);
		}
	}
	FORI(n) {
		int t[2] = {p[i][0], p[i][1]};
		int z[2] = {p[i][2], p[i][3]};
		int Z[2] = {p[i][4], p[i][5]};
		if (inner_point(s[0], s[1], t[0], t[1]) &&
				inner_point(s[0], s[1], z[0], z[1]) &&
				inner_point(s[0], s[1], Z[0], Z[1]))
		{
			o[t[1]][t[0]] = RGBA_YELLOW;
			o[z[1]][z[0]] = RGBA_GREEN;
			o[Z[1]][Z[0]] = RGBA_MAGENTA;
		}
	}
	iio_save_image_uint8_vec("-", (uint8_t*)o, s[0], s[1], 4);
	return EXIT_SUCCESS;
}

static void put_pixel(int a, int b, void *pp)
{
	struct {int w, h; struct rgba_value *x, c;} *p = pp;
	if (inner_point(p->w, p->h, a, b))
		p->x[b*p->w+a] = p->c;
}

static void overlay_segment_color(struct rgba_value *x, int w, int h,
		int px, int py, int qx, int qy, struct rgba_value c)
{
	struct {int w, h; struct rgba_value *x, c; } e = {w, h, x, c};
	traverse_segment(px, py, qx, qy, put_pixel, &e);
}

static void overlay_line(double a, double b, double c,
		struct rgba_value *x, int w, int h, struct rgba_value k)
{
	if (b == 0) {
		int f[2] = {-c/a, 0};
		int t[2] = {-c/a, h-1};
		overlay_segment_color(x, w, h, f[0], f[1], t[0], t[1], k);
	} else {
		double alpha = -a/b;
		double beta = -c/b;
		int f[2] = {0, beta};
		int t[2] = {w-1, alpha*(w-1)+beta};
		overlay_segment_color(x, w, h, f[0], f[1], t[0], t[1], k);
	}
}

static void epipolar_line(double L[3], double F[9], double x[2])
{
	L[0] = x[0]*F[0] + x[1]*F[3] + F[6];
	L[1] = x[0]*F[1] + x[1]*F[4] + F[7];
	L[2] = x[0]*F[2] + x[1]*F[5] + F[8];
}

static void project_point_to_line(double px[2], double L[3], double x[2])
{
	double a = L[0];
	double b = L[1];
	double c = L[2];
	double d = a*a + b*b;
	double e = a*x[1] - b*x[0];
	px[0] = (-a*c - b*e)/d;
	px[1] = (-b*c + a*e)/d;
}

int main_viewepi(int c, char *v[])
{
	if (c != 12 && c != 13) {
		fprintf(stderr, "usage:\n\t%s f0 ... f8 sx sy [mask]\n", *v);
		//                          0 1      9  10 11  12
		return EXIT_FAILURE;
	}
	double A[9]; FORI(9) A[i] = atof(v[1+i]);
	int s[2], n; FORI(2) s[i] = atoi(v[10+i]);
	double (*p)[4] = (void*)read_ascii_doubles(stdin, &n);
	n /= 4;
	struct rgba_value (*o)[s[0]] = xmalloc(s[0]*s[1]*4);
	struct {int w, h; struct rgba_value *x; } usr = {s[0], s[1], o[0]};
	bool bmask[n]; FORI(n) bmask[i] = c<13;
	if (c>12) { FILE *f = xfopen(v[12], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=t;}
		xfclose(f);
	}
	FORI(s[0]*s[1]) o[0][i] = RGBA_BLACK;
	// for each point p (=px) in image A
	// 1. draw the epipolar line L of p
	FORI(n) { double L[3]; epipolar_line(L, A, p[i]);
		if (!bmask[i])
			overlay_line(L[0],L[1],L[2],o[0],s[0],s[1],RGBA_GRAY10);
	}
	FORI(n) { double L[3]; epipolar_line(L, A, p[i]);
		if (bmask[i])
			overlay_line(L[0],L[1],L[2],o[0],s[0],s[1],RGBA_GRAY50);
	}
	// 3. draw the projection line of q to L
	FORI(n) {
		double L[3];
		epipolar_line(L, A, p[i]);
		double Ly[2];
		int y[2] = {lrint(p[i][2]), lrint(p[i][3])};
		project_point_to_line(Ly, L, p[i]+2);
		int iLy[2] = {lrint(Ly[0]), lrint(Ly[1])};
		if (inner_point(s[0], s[1], iLy[0], iLy[1]))
		{
			if (inner_point(s[0], s[1], y[0], y[1]))
				traverse_segment(y[0], y[1], iLy[0], iLy[1],
						bmask[i]?draw_brighter_pixel:
						draw_phantom_pixel, &usr);
			if (bmask[i])
				o[iLy[1]][iLy[0]] = RGBA_MAGENTA;
			else
				o[iLy[1]][iLy[0]] = RGBA_RED;
		}
	}
	// 2. draw the corresponding point q (py)
	FORI(n) {
		int y[2] = {lrint(p[i][2]), lrint(p[i][3])};
		if (inner_point(s[0], s[1], y[0], y[1])) {
			if (bmask[i])
				o[y[1]][y[0]] = RGBA_GREEN;
			else
				o[y[1]][y[0]] = RGBA_BLUE;
		}
	}
	iio_save_image_uint8_vec("-", (uint8_t*)o, s[0], s[1], 4);
	return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
	assert(4 == sizeof(struct rgba_value));
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "points")) return main_viewp(c-1, v+1);
	else if (0 == strcmp(v[1], "pairs")) return main_viewpairs(c-1, v+1);
	else if (0 == strcmp(v[1], "triplets")) return main_viewtrips(c-1, v+1);
	else if (0 == strcmp(v[1], "epipolar")) return main_viewepi(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s [points|pairs|triplets|epipolar] "
			       "params... < data.txt | display\n", *v);
	       return EXIT_FAILURE;
	}
}
