// block matching of two images guided by a fundamental matrix
// (without rectification)

// TODO:
// 1. optimize getpixel by omitting border computations
// 2. inline LRRL
// 3. multi-scale
// 4. inline flath
// 5. inline mindiff

#include <assert.h>
#include <math.h>


static int window_image_square_5x5[] = {5,5, 2,2,
	1,1,1,1,1,
	1,1,1,1,1,
	1,1,2,1,1,
	1,1,1,1,1,
	1,1,1,1,1,
};
static int window_image_horiz_25[] = {9,3, 4,1,
	1,1,1,1,0,1,1,1,1,
	1,1,1,1,2,1,1,1,1,
	1,1,1,1,0,1,1,1,1,
};
static int window_image_vert_25[] = {3,9, 1,4,
	1,1,1,
	1,1,1,
	1,1,1,
	1,1,1,
	0,2,0,
	1,1,1,
	1,1,1,
	1,1,1,
	1,1,1,
};
static int window_image_diag_25[] = {9,9, 4,4,
	0,0,0,0,0,0,0,1,1,
	0,0,0,0,0,0,1,1,1,
	0,0,0,0,0,1,1,1,0,
	0,0,0,0,1,1,1,0,0,
	0,0,0,1,2,1,0,0,0,
	0,0,1,1,1,0,0,0,0,
	0,1,1,1,0,0,0,0,0,
	1,1,1,0,0,0,0,0,0,
	1,1,0,0,0,0,0,0,0,
};
static int window_image_ndiag_25[] = {9,9, 4,4,
	1,1,0,0,0,0,0,0,0,
	1,1,1,0,0,0,0,0,0,
	0,1,1,1,0,0,0,0,0,
	0,0,1,1,1,0,0,0,0,
	0,0,0,1,2,1,0,0,0,
	0,0,0,0,1,1,1,0,0,
	0,0,0,0,0,1,1,1,0,
	0,0,0,0,0,0,1,1,1,
	0,0,0,0,0,0,0,1,1,
};
static int *window_images[] = {
		window_image_square_5x5,
		window_image_horiz_25,
		window_image_vert_25,
		window_image_diag_25,
		window_image_ndiag_25,
};
static int number_of_window_images = sizeof(window_images)/sizeof*window_images;


#include "xmalloc.c"
#include "getpixel.c"

static int plot_centered_segment(int (*P)[2], int w, int h,
		double a, double b, double c,
		double px, double py, double rad
		)
{
	if (!isfinite(rad)) rad = fmax(w,h);
	assert(fabs(a*px + b*py + c) < 10e-6);
	double dirn = rad/hypot(a,b), sgn = a - b > 0 ? 1 : -1;
	double dir[2] = {-b*dirn*sgn, a*dirn*sgn};
	assert(dir[0] + dir[1] > 0);
	int r = 0;
	if (fabs(a) < fabs(b)) { // slope less than 1
		double p = -a/b;
		double q = -c/b;
		int xmin = round(fmax(0,   px - dir[0]));
		int xmax = round(fmax(w-1, px + dir[0]));
		if (xmin > xmax) return 0;
		for (int x = xmin; x <= xmax; x++) {
			int y = round(p*x + q);
			P[r][0] = x;
			if (y >= 0 && y < h)
				P[r++][1] = y;
		}
	} else {
		double p = -b/a;
		double q = -c/a;
		int ymin = round(fmax(0,   py - dir[1]));
		int ymax = round(fmin(h-1, py + dir[1]));
		if (ymin > ymax) return 0;
		for (int y = ymin; y <= ymax; y++) {
			int x = round(p*y + q);
			P[r][1] = y;
			if (x >= 0 && x < w)
				P[r++][0] = x;
		}
	}
	return r;
}

static void project_point_to_line(double pp[2],
		double a, double b, double c, double x, double y)
{
	pp[0] = (b*b*x - a*b*y - a*c)/(a*a + b*b);
	pp[1] = (a*a*y - a*b*x - b*c)/(a*a + b*b);
}

// Observation: it is unnecessary to optimize this function.
// According to callgrind, the program spends less than 0.1%
// of the overall running time inside it.
static int plot_epipolar_fancy(int (*p)[2], double fm[9], int w, int h,
		int i, int j, float ini[2], float rad)
{
	double a = i*fm[0] + j*fm[3] + fm[6];
	double b = i*fm[1] + j*fm[4] + fm[7];
	double c = i*fm[2] + j*fm[5] + fm[8];
	double ij[2] = {i, j};
	if (isfinite(*ini)) {
		ij[0] += ini[0];
		ij[1] += ini[1];
	}
	double pp[2];
	project_point_to_line(pp, a, b, c, ij[0], ij[1]);
	return plot_centered_segment(p, w, h, a, b, c, pp[0], pp[1], rad);
}

typedef float (*vector_correlation_measure)(float*,float*,int);

// Observation: it is unnecessary to optimize this function.
// According to callgrind, the program spends about 5% of the running time
// inside it
static float ssd_minus_mean(float *x, float *y, int n)
{
	float mx = 0, my = 0;
	for (int i = 0; i < n; i++)
	{
		mx += x[i]/n;
		my += y[i]/n;
	}

	float r = 0;
	for (int i = 0; i < n; i++)
	{
		float s = (x[i] - mx) - (y[i] - my);
		r += s * s;
	}
	return r;
}

struct correlation_window {
	int n;
	int (*off)[2];
};

struct correlation_window *create_window_list(int *nwin)
{
	int n = number_of_window_images;
	int **t = window_images;
	struct correlation_window *r = xmalloc((1+n)*sizeof*r);
	for (int i = 0; i < n; i++) {
		struct correlation_window *win = r + i;
		int w = t[i][0];
		int h = t[i][1];
		int cx = t[i][2];
		int cy = t[i][3];
		int count = 0;
		for (int p = 0; p < w*h; p++)
			if (t[i][4+p])
				count += 1;
		win->n = count;
		win->off = xmalloc(count*2*sizeof(int));
		for (int q = 0; q < h; q++)
		for (int p = 0; p < w; p++)
			if (t[i][4+w*q+p]) {
				count -= 1;
				win->off[count][0] = p - cx;
				win->off[count][1] = q - cy;
				if (2 == t[i][4+w*q+p])
					assert(p == cx && q == cy);
			}
		assert(count == 0);
		assert(2 == t[i][4+w*cy+cx]);
	}
	r[n].n = 0;
	*nwin = n;
	return r;
}

static void print_window_list(struct correlation_window *t)
{
	for (int i = 0; ; i++)
		if (t[i].n) {
			fprintf(stderr, "window %d on list has %d offsets\n",
					i, t[i].n);
			for (int j = 0 ; j < t[i].n; j++)
				fprintf(stderr, "\toff[%d] = %d %d\n",
						j, t[i].off[j][0], t[i].off[j][1]);
			fprintf(stderr, "\n");
		}
		else
			break;
}

static void free_window_list(struct correlation_window *t)
{
	for (int i = 0; ; i++)
		if (t[i].n)
			free(t[i].off);
		else
			break;
	free(t);
}

static struct correlation_window *global_table_of_windows = NULL;
static int global_number_of_windows = 0;


// this function is where 95% of the running time is spent
// (inside the loop)
static double corr(float *a, float *b, int w, int h, int pd,
		int ax, int ay, int bx, int by, int wintype)
{
	struct correlation_window *wen = global_table_of_windows+wintype;
	int n = pd * wen->n;
	float pa[n], pb[n];
	int cx = 0;
	for (int i = 0; i < wen->n; i++)
	for (int l = 0; l < pd; l++)
	{
		int dx = wen->off[i][0];
		int dy = wen->off[i][1];
		pa[cx] = getsample_2(a, w, h, pd, ax + dx, ay + dy, l);
		pb[cx] = getsample_2(b, w, h, pd, bx + dx, by + dy, l);
		cx += 1;
	}

	vector_correlation_measure f = ssd_minus_mean;
	return f(pa, pb, n);
}

void bmfm_fancy(float *disp,         // output disparities image (dx, dy)
		float *errc,         // output error image
		float *a,            // input image A
		float *b,            // input image B
		int w,               // width
		int h,               // height
		int pd,              // pixel dimension
		double fm[9],        // fundamental matrix
		float *disp_init,    // initialization (optional)
		float *search_radius // optional, w.r.t. initialization
		)
{
	int maxpoints = 2 * (w+h), (*P)[2] = xmalloc(maxpoints*sizeof*P);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float rad = NAN, ini[2] = {0, 0};
		if (search_radius) rad = search_radius[j*w+i];
		if (disp_init) ini[0] = disp_init[2*(j*w+i)+0];
		if (disp_init) ini[1] = disp_init[2*(j*w+i)+1];
		int np = plot_epipolar_fancy(P, fm, w, h, i, j, ini, rad);
		float mincorr = INFINITY;
		int minidx = 0;
		for (int k = 0; k < np; k++)
		for (int l = 0; l < global_number_of_windows; l++)
		{
			int p = P[k][0];
			int q = P[k][1];
			float c = corr(a, b, w, h, pd, i, j, p, q, l);
			if (c < mincorr) {
				mincorr = c;
				minidx = k;
			}
		}
		if (disp) disp[2*(j*w+i) + 0] = P[minidx][0] - i;
		if (disp) disp[2*(j*w+i) + 1] = P[minidx][1] - j;
		if (errc) errc[j*w+i] = mincorr;
	}
	free(P);
}

#ifdef MAIN_BMFM

#include <stdio.h>
#include <string.h>
#include "iio.h"

#include "fail.c"
#include "parsenumbers.c"

int main(int c, char *v[])
{
	if (c != 8) {
		fprintf(stderr, "usage:\n\t"
			"%s a.png b.png \"fm\" out_disp out_err ini rad\n", *v);
		//        0 1     2       3    4        5       6   7
		return 1;
	}
	char *filename_a = v[1];
	char *filename_b = v[2];
	char *fm_text = v[3];
	char *filename_out_disp = v[4];
	char *filename_out_err = v[5];
	char *filename_init = v[6];
	float maxradius = atof(v[7]);

	int nfm;
	float *ffm = alloc_parse_floats(9, fm_text, &nfm);
	if (nfm != 9)
		fail("expects a fundamental matrix (9 numbers)");
	double fm[9];
	for (int i = 0; i < 9; i++)
		fm[i] = ffm[i];
	free(ffm);

	int w, h, pd, ww, hh, ppdd;
	float *a = iio_read_image_float_vec(filename_a, &w, &h, &pd);
	float *b = iio_read_image_float_vec(filename_b, &ww, &hh, &ppdd);
	if (w != ww || h != hh || pd != ppdd)
		fail("input images size mismatch");

	float *o = xmalloc(w*h*3*sizeof*o);
	float *e = xmalloc(w*h*sizeof*o);
	float *rad = xmalloc(w*h*sizeof*o);
	for (int i = 0; i < w*h; i++)
		rad[i] = maxradius;

	float *i = NULL;
	if (0 != strcmp(filename_init, "0")) {
		i = iio_read_image_float_vec(filename_init, &ww, &hh, &ppdd);
		if (w != ww || h != hh)
			fail("init image size mismatch");
		if (ppdd != 2)
			fail("init file must be a vector field");
	}

	global_table_of_windows = create_window_list(&global_number_of_windows);
	//print_window_list(global_table_of_windows);
	global_number_of_windows = 1;

	bmfm_fancy(o, e, a, b, w, h, pd, fm, i, rad);

	iio_save_image_float_vec(filename_out_disp, o, w, h, 2);
	iio_save_image_float(filename_out_err, e, w, h);

	free(a);
	free(b);
	free(o);
	free(e);
	if (i) free(i);

	return 0;
}

#endif//MAIN_BMFM
