#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"
#include "marching_squares.c"
#include "marching_interpolation.c"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

// print an error message and abort the program
static void error(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");

#ifdef NDEBUG
	exit(-1);
#else
	exit(*(int *)0x43);
#endif//NDEBUG
}

//// always returns a valid pointer
//static void *xmalloc(size_t size)
//{
//	if (size == 0)
//		error("xmalloc: zero size");
//	void *new = malloc(size);
//	if (!new)
//	{
//		double sm = size / (0x100000 * 1.0);
//		error("xmalloc: out of memory when requesting "
//				"%zu bytes (%gMB)",//:\"%s\"",
//				size, sm);//, strerror(errno));
//	}
//	return new;
//}


//// utility function: alloc a 2d matrix contiguously (wxh elements of size n)
//static void *matrix_build(int w, int h, size_t n)
//{
//	size_t p = sizeof(void *);
//	char *r = xmalloc(h*p + w*h*n);
//	for (int i = 0; i < h; i++)
//		*(void **)(r + i*p) = r + h*p + i*w*n;
//	return r;
//}

// visualization parameters that can be applied to an arbitrary image
struct closeup_params {
	float black;
	float white;
	int cropbox[4];
	int zoom;
	int zoomtype;
	int number_of_level_lines;
	float *levels;
	char *level_type;
};

static int bound(int a, int x, int b)
{
	fprintf(stderr, "bound %d %d %d\n", a, x, b);
	assert(a <= b);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

//static float fbound(float a, float x, float b)
//{
//	assert(a <= b);
//	if (x < a) return a;
//	if (x > b) return b;
//	return x;
//}

static float *crop(float *x, int w, int h, int pd, int c[4], int *ow, int *oh)
{
	int x0 = bound(0, c[0], w-1);
	int y0 = bound(0, c[1], h-1);
	int xf = bound(x0+1, c[2]==-1?w:c[2], w);
	int yf = bound(y0+1, c[3]==-1?h:c[3], h);
	int cw = xf - x0;
	int ch = yf - y0;
	float *y = xmalloc(cw*ch*pd*sizeof*y);
	FORJ(ch) FORI(cw) FORL(pd)
		y[(i + j*cw)*pd + l] = x[(x0 + i + (y0 + j)*w)*pd + l];
	*ow = cw;
	*oh = ch;
	return y;
}

//static float (**crop(float (**x)[3], int w, int h, int c[4],
//		int *ow, int *oh))[3]
//{
//	int x0 = bound(0, c[0], w-1);
//	int y0 = bound(0, c[1], h-1);
//	int xf = bound(x0+1, c[2], w);
//	int yf = bound(y0+1, c[3], h);
//	int cw = xf - x0;
//	int ch = yf - y0;
//	float (**y)[3] = matrix_build(cw, ch, sizeof**y);
//	FORJ(ch) FORI(cw) FORL(3)
//		y[j][i][l] = x[j+y0][i+x0][l];
//	*ow = cw;
//	*oh = ch;
//	return y;
//}

#define BAD_MIN(a,b) (b)<(a)?(b):(a)

static void assert_rgb(double t[3])
{
	for (int i = 0; i < 3; i++)
		assert(t[i] >= 0 && t[i] <= 1);
}

static void assert_hsv(double t[3])
{
	if (t[0] < 0 || t[0] >= 360) error("queca %g\n", t[0]);
	assert(t[0] >= 0 && t[0] < 360);
	if (!(t[1] >= 0 && t[1] <= 1))
		error("CACA S = %g\n", t[1]);
	assert(t[1] >= 0 && t[1] <= 1);
	assert(t[2] >= 0 && t[2] <= 1);
}

static void rgb_to_hsv_doubles(double *out, double *in)
{
	assert_rgb(in);
	double r, g, b, h, s, v, M, m;
	r = in[0]; g = in[1]; b = in[2];

	//printf("rgb %g,%g,%g...\n", r, g, b);

	if (g >= r && g >= b) {
		M = g;
		m = BAD_MIN(r, b);
		h = M == m ? 0 : 60*(b-r)/(M-m)+120;
	}
	else if (b >= g && b >= r) {
		M = b;
		m = BAD_MIN(r, b);
		h = M == m ? 0 : 60*(r-g)/(M-m)+240;
	}
	else {
		assert(r >= g && r >= b);
		M = r;
		if (g >= b) {
			m = b;
			h = M == m ? 0 : 60*(g-b)/(M-m)+0;
		} else {
			m = g;
			h = M == m ? 0 : 60*(g-b)/(M-m)+360;
		}
	}

	s = M == 0 ? 0 : (M - m)/M;
	v = M;
	h = fmod(h, 360);

	//printf("\thsv %g,%g,%g\n", h, s, v);
	out[0] = h; out[1] = s; out[2] = v;
	assert_hsv(out);
}

static void get_hsv(float *outh, float *outs, float *outv, float x[3])
{
	double din[3] = {x[0]/255.0, x[1]/255.0, x[2]/255.0};
	double dout[3];
	rgb_to_hsv_doubles(dout, din);
	if (outh) *outh = 255.0*dout[0]/360.0;
	if (outs) *outs = 255.0*dout[1];
	if (outv) *outv = 255.0*dout[2];
}

static float get_rgb_intensity(float x[3])
{
	return 0.299 * x[0] + 0.587 * x[1] + 0.114 * x[2];
}

// TODO: add more scalarization options here: CIELABs, XYZ, PCA ...
static float *scalar(float *x, int w, int h, int pd, char *level_type)
{
	if (pd == 1) return x;
	float *y = xmalloc(w * h * pd * sizeof*y);
	if (false) { ;
	} else if (pd == 3 && 0 == strcmp(level_type, "h")) {
		FORI(w*h) get_hsv(y+i, NULL, NULL, x + i*pd);
	} else if (pd == 3 && 0 == strcmp(level_type, "s")) {
		FORI(w*h) get_hsv(NULL, y+i, NULL, x + i*pd);
	} else if (pd == 3 && 0 == strcmp(level_type, "v")) {
		FORI(w*h) get_hsv(NULL, NULL, y+i, x + i*pd);
	} else if (pd == 3 && 0 == strcmp(level_type, "i")) {
		FORI(w*h) y[i] = get_rgb_intensity(x + i*pd);
	} else if (pd == 3 && 0 == strcmp(level_type, "r")) {
		FORI(w*h) y[i] = x[i*pd];
	} else if (pd == 3 && 0 == strcmp(level_type, "g")) {
		FORI(w*h) y[i] = x[i*pd + 1];
	} else if (pd == 3 && 0 == strcmp(level_type, "b")) {
		FORI(w*h) y[i] = x[i*pd + 2];
	} else if (isdigit(level_type[0])) {
		int c = bound(0, atoi(level_type), pd-1);
		FORI(w*h) y[i] = x[i*pd + c];
	} else {
		FORI(w*h) {
			float m = 0;
			FORL(pd) m += x[i*pd + l];
			y[i] = m/pd;
		}
	}
	return y;
}

static float getsample(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

static void setsample(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

static float interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	// return a;
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return interpolate_nearest(a, b, c, d, x, y);
	case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	default: error("caca de vaca");
	}
	return -1;
}

static void interpolate_vec(float *out, float *x, int w, int h, int pd,
		float p, float q, int m)
{
	if (p < 0 || q < 0 || p+1 >= w || q+1 >= h) {
		FORL(pd) out[l] = 0;
	} else {
		int ip = floor(p);
		int iq = floor(q);
		FORL(pd) {
			float a = getsample(x, w, h, pd, ip  , iq  , l);
			float b = getsample(x, w, h, pd, ip  , iq+1, l);
			float c = getsample(x, w, h, pd, ip+1, iq  , l);
			float d = getsample(x, w, h, pd, ip+1, iq+1, l);
			float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
			//fprintf(stderr, "p%g q%g ip%d iq%d a%g b%g c%g d%g l%d v%g\n", p, q, ip, iq, a, b, c, d, l, v);
			out[l] = v;
		}
	}
}

static float *zoom(float *x, int w, int h, int pd, int n, int zt,
		int *ow, int *oh)
{
	int W = n*w - n + 1;
	int H = n*h - n + 1;
	float *y = xmalloc(W*H*pd*sizeof*y), nf = n;
	FORJ(H) FORI(W) {
		float tmp[pd];
		interpolate_vec(tmp, x, w, h, pd, i/nf, j/nf, zt);
		FORL(pd)
			setsample(y, W, H, pd, i, j, l, tmp[l]);
	}

	*ow = W;
	*oh = H;
	return y;
}

// draw a segment between two points
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px < qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py < qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px+j*slope), j+py, e);
		}
	}
}

static int hack_width;
static int hack_height;

static void red_fpixel(int x, int y, void *ii)
{
	if (x < 0 || y < 0 || x > hack_width || y > hack_height) return;
	//float (**i)[3] = ii;
	//i[y][x][0] = 255;
	//i[y][x][1] = 0;
	//i[y][x][2] = 0;
	float *i = ii;
	setsample(i, hack_width, hack_height, 3, x, y, 0, 255);
	setsample(i, hack_width, hack_height, 3, x, y, 1, 0);
	setsample(i, hack_width, hack_height, 3, x, y, 2, 0);
}

static void overlay_level_line_in_red(float *y, int n,
		float *x, int w, int h, float t)
{
	hack_width = n*(w-1)+1;
	hack_height = n*(h-1)+1;
	int ns;
	float (*s)[2][2] = marching_squares_whole_image_float(&ns, x, w, h, t);
	for (int i = 0; i < ns; i++) {
		traverse_segment(n*s[i][0][0], n*s[i][0][1],
				n*s[i][1][0], n*s[i][1][1],
				red_fpixel, y);
	}
	free(s);
}

static float *closeup(float *x, int w, int h, int pd, struct closeup_params *p,
		int *ow, int *oh, int *opd)
{
	// crop the image
	int cw, ch;
	float *cx = crop(x, w, h, pd, p->cropbox, &cw, &ch);
	iio_save_image_float_vec("/tmp/cropped.png", cx, cw, ch, pd);

	// compute the scalar channel
	float *scx = scalar(cx, cw, ch, pd, p->level_type);
	iio_save_image_float_vec("/tmp/scalar.png", scx, cw, ch, 1);

	// compute the zoomed image
	int zcw, zch;
	float *zcx = cx;
	//if (p->zoom > 1)
		zcx = zoom(cx, cw, ch, pd, p->zoom, p->zoomtype, &zcw, &zch);

	// create output image
	int out_w = zcw, out_h = zch, out_pd = pd;
	float *out_x = zcx;
	if (p->number_of_level_lines > 0 && pd == 3)
		out_x = zcx;
	if (p->number_of_level_lines > 0 && pd != 3) {
		out_pd = 3;
		out_x = xmalloc(out_w * out_h * 3 * sizeof*out_x);
		FORJ(out_h) FORI(out_w) FORL(out_pd) {
			float g = l<pd ?
				getsample(zcx,zcw,zch,pd, i,j,l) :
				getsample(zcx,zcw,zch,pd, i,j,pd-1);
			setsample(out_x,out_w,out_h,out_pd, i,j,l, g);
		}
	}

	// adjust contrast of output image
	if (p->white || p->black) {
		FORJ(out_h) FORI(out_w) FORL(out_pd) {
			float g = getsample(out_x,out_w,out_h,out_pd, i,j,l);
			g = floor(255 * (g - p->black)/(p->white - p->black));
			if (g < 0) g = 0;
			if (g > 255) g = 255;
			setsample(out_x,out_w,out_h,out_pd, i,j,l, g);
		}
	}

	// overlay level lines of "scalar" on output image
	//if (p->number_of_level_lines) {
	//}
	FORI(p->number_of_level_lines) {
		assert(out_pd == 3);
		//fprintf(stderr, "pzoom = %d\n", p->zoom);
		//fprintf(stderr, "cw = %d\n", cw);
		//fprintf(stderr, "out_w = %d\n", out_w);
		assert(p->zoom*(cw - 1) + 1 == out_w);
		assert(p->zoom*(ch - 1) + 1 == out_h);
		float lev = p->levels[i];
		overlay_level_line_in_red(out_x, p->zoom, scx, cw, ch, lev);
	}

	// free temporary stuff

	*ow = out_w;
	*oh = out_h;
	*opd = out_pd;

	return out_x;
}

//static int parse_integers(int *t, int nmax, const char *s)
//{
//	int i = 0, w;
//	while (i < nmax && 1 == sscanf(s, "%d %n", t + i, &w)) {
//		i += 1;
//		s += w;
//	}
//	return i;
//}

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

int main(int c, char *v[])
{
	if (c != 13) {
		fprintf(stderr, "usage:\n\t%s "
	"black white cropx0 cropy0 cropw croph zoom zoomtype \"levels\" levchan"
	// 1    2      3      4      5    6      7   8         9         10
	" input_image output_image" "\n" , *v);
	// 11              12
		return EXIT_FAILURE;
	}

	struct closeup_params p[1];

	p->black = atof(v[1]);
	p->white = atof(v[2]);
	FORI(4) p->cropbox[i] = atof(v[3+i]);
	p->zoom = atoi(v[7]);
	p->zoomtype = atoi(v[8]);
	//p->number_of_level_lines = 0; // parse v[9]
	//p->levels = NULL;             // parse v[9]
	float pi[100]; p->levels = pi;
	p->number_of_level_lines = parse_floats(pi, 100, v[9]);
	//int rpi = parse_floats(pi, 10, v[9]);
	//fprintf(stderr, "parsed %d floats from \"%s\":\n", rpi, v[9]);
	//FORI(rpi) fprintf(stderr, "float[%d] = %g\n", i, pi[i]);
	p->level_type = v[10];

	int w, h, pd;
	float *x = iio_read_image_float_vec(v[11], &w, &h, &pd);

	int ow, od, opd;
	float *y = closeup(x, w, h, pd, p, &ow, &od, &opd);

	iio_save_image_float_vec(v[12], y, ow, od, opd);

	free(x);
	free(y);

	return EXIT_SUCCESS;
}
