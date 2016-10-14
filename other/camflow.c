/**
 * Display video from webcam
 *
 * Author  Nash
 * License GPL
 * Website http://nashruddin.com
 */

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "cv.h"
#include "highgui.h"



#include "iio.h"
#include "harris/harris.c"
#include "harris/gauss.c"
#include "harris/image.c"
#include "harris/ntuple.c"
#include "harris/misc.c"

#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfont9x15.c"
#include "fonts/xfont8x13.c"

#include "seconds.c"


char *global_method_id;
float global_scale_step;
int   global_nscales;
int   global_last_scale;
float global_params[0x100];
int   global_nparams;
float global_vscale;
const float global_vscale_factor = 1.4142;

struct bitmap_font global_font;

double global_harris_sigma = 1;    // s
double global_harris_k = 0.04;     // k
double global_harris_flat_th = 20; // t
int    global_harris_neigh = 3;    // n

enum flow_visualization_id {
	VFLOW_COLORS,
	VFLOW_ARROWS,
	VFLOW_DIVERGENCE,
	VFLOW_GRADIENT,
	VFLOW_BACK,
	VFLOW_BACKDIFF,
	VFLOW_BARROWS,
};

int global_flow_visualization = VFLOW_COLORS;

bool global_display_img = false;
bool global_display_diff = false;

#include "xmalloc.c"

static void hsv_to_rgb_doubles(double *out, double *in)
{
	//assert_hsv(in);
	double r, g, b, h, s, v; r=g=b=h=s=v=0;
	h = in[0]; s = in[1]; v = in[2];
	if (s == 0)
		r = g = b = v;
	else {
		int H = fmod(floor(h/60),6);
		double p, q, t, f = h/60 - H;
		p = v * (1 - s);
		q = v * (1 - f*s);
		t = v * (1 - (1 - f)*s);
		switch (H) {
			case 6:
			case 0: r = v; g = t; b = p; break;
			case 1: r = q; g = v; b = p; break;
			case 2: r = p; g = v; b = t; break;
			case 3: r = p; g = q; b = v; break;
			case 4: r = t; g = p; b = v; break;
			case -1:
			case 5: r = v; g = p; b = q; break;
			default:
				fprintf(stderr, "H=%d\n", H);
				assert(false);
		}
	}
	out[0] = r; out[1] = g; out[2] = b;
	//assert_rgb(out);
}

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

static void flow_to_frgb(float *y, float u, float v, float m)
{
	double r = hypot(u, v);
	if (r > 1e3 || !isfinite(r))
		y[0] = y[1] = y[2] = 255;
	r = r>m ? 1 : r/m;
	double a = atan2(v, -u);
	a = (a+M_PI)*(180/M_PI);
	a = fmod(a, 360);
	double hsv[3], rgb[3];
	hsv[0] = a;
	hsv[1] = r;
	hsv[2] = r;
	hsv_to_rgb_doubles(rgb, hsv);
	for (int l = 0; l < 3; l++)
		y[l] = 255*rgb[l];
}

static void flow_to_farr(float *farr, float *u, float *v, int w, int h, float m)
{
	float *f = xmalloc(2*w*h*sizeof*f);
	float *a = xmalloc(w*h*sizeof*a);
	for (int i = 0; i < w*h; i++) {
		f[2*i+0] = u[i];
		f[2*i+1] = v[i];
		a[i] = 255;
	}
	void flowarrows(float *, float *, int, int, float, int);
	flowarrows(a, f, w, h, 0.037*m, 19);
	for (int i = 0; i < w*h; i++)
		for (int l = 0; l < 3; l++)
			farr[3*i+l] = a[i];
	free(f);
	free(a);
}

static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

static void flow_to_fdiv(float *fdiv, float *u, float *v, int w, int h, float m)
{
	float (*cdiv)[w][3] = (void*)fdiv;
#define FX(i,j) getpixel(u,w,h,(i),(j))
#define FY(i,j) getpixel(v,w,h,(i),(j))
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float xx = -FX(i-1,j-1)-2*FX(i-1,j)-FX(i-1,j+1)
			+FX(i+1,j-1)+2*FX(i+1,j)+FX(i+1,j+1);
		float yy = -FY(i-1,j-1)-2*FY(i,j-1)-FY(i+1,j-1)
			+FY(i-1,j+1)+2*FY(i,j+1)+FY(i+1,j+1);
		float dd = xx + yy;
		float g = 255.0*(dd + m)/(2.0*m);
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		for (int l = 0; l < 3; l++)
			cdiv[j][i][l] = g;
	}
#undef FX
#undef FY
}

static float getsamplen(float *fx, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	float (*x)[w][pd] = (void*)fx;
	return x[j][i][l];
	//return x[(i+j*w)*pd + l];
}

static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}


static void bilinear_interpolation_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q)
{
	int ip = p;
	int iq = q;
	for (int l = 0; l < pd; l++) {
		float a = getsamplen(x, w, h, pd, ip  , iq  , l);
		float b = getsamplen(x, w, h, pd, ip+1, iq  , l);
		float c = getsamplen(x, w, h, pd, ip  , iq+1, l);
		float d = getsamplen(x, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}

static void flow_to_fback(float *wb, float *b, float *u, float *v, int w, int h)
{
	float (*out)[w][3] = (void*)wb;
	//float (*in)[w][3] = (void*)b;
	float (*img_u)[w] = (void*)u;
	float (*img_v)[w] = (void*)v;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		float p[2] = {i + img_u[j][i], j + img_v[j][i]};
		float result[3];
		bilinear_interpolation_at(result, b, w, h, 3, p[0], p[1]);
		float factor = 1;
		for (int l = 0; l < 3; l++)
			out[j][i][l] = result[l];
	}
}

static void flow_to_fbackdiff(float *out, float *a, float *b,
		float *u, float *v, int w, int h, float m)
{
	float *wb = xmalloc(w * h * 3 * sizeof*wb);
	flow_to_fback(wb, b, u, v, w, h);
	for (int i = 0; i < w * h * 3; i++) {
		float dd = wb[i] - a[i];
		float g = dd/1 + 128;
		if (g < 0) g = 0;
		if (g > 255) g = 255;
		out[i] = g;
	}
	free(wb);
}

static void flow_to_barr(float *barr, float *a, float *u, float *v,
		int w, int h, float m)
{
	float *f = xmalloc(2*w*h*sizeof*f);
	float *arr = xmalloc(w*h*sizeof*a);
	for (int i = 0; i < w*h; i++) {
		f[2*i+0] = u[i];
		f[2*i+1] = v[i];
		arr[i] = a[3*i+1]/2+127;
	}
	void flowarrows(float *, float *, int, int, float, int);
	flowarrows(arr, f, w, h, 0.037*m, 19);
	for (int i = 0; i < w*h; i++)
		if (arr[i] > 127)
			for (int l = 0; l < 3; l++)
				barr[3*i+l] = a[3*i+l]/2+127;
		else
			for (int l = 0; l < 3; l++)
				barr[3*i+l] = arr[i];
	free(f);
	free(arr);
}

void multi_scale_optical_flow_pd(char *algorithm_name, float *pars, int npars,
		float *u, float *v, float *a, float *b, int w, int h, int pd,
		int nscales, float scale_step, int last_scale);

// visualize the "difference" between two consecutive frames
static void d_tacu(float *out, float *in_a, float *in_b, int w, int h, int pd)
{
	if (global_display_img) {
		for (int i = 0; i < w*h*pd; i++) {
			float g = in_b[i];
			if (g > 255) g = 255;
			if (g < 0) g = 0;
			out[i] = g;
		}
		return;
	}

	if (global_display_diff) {
		for (int l = 0; l < pd; l++)
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			float a00 = in_a[(j*w+i)*pd+l];
			float a10 = in_a[(j*w+i+1)*pd+l];
			float a01 = in_a[(j*w+i+w)*pd+l];
			float adx = a10 - a00;
			float ady = a01 - a00;
			float na2 = adx*adx + ady*ady;
			float b00 = in_b[(j*w+i)*pd+l];
			float idt = b00 - a00;
			float g = 127 + idt;//*na2/100;
			//if (idt > 25) g = 255;
			//if (idt < -25) g = 0;
			if (g > 255) g = 255;
			if (g < 0) g = 0;
			out[(j*w+i)*pd+l] = g;
		}
		//for (int i = 0; i < w*h*pd; i++) {
		//	float g = (in_b[i] - in_a[i])/1 + 128;
		//	if (g > 255) g = 255;
		//	if (g < 0) g = 0;
		//	out[i] = g;
		//}
		return;
	}


	float *u = xmalloc(w*h*sizeof*u);
	float *v = xmalloc(w*h*sizeof*u);

	if (0 == strcmp(global_method_id, "tvl1")) {
		float *I0 = xmalloc(w*h*sizeof*I0);
		float *I1 = xmalloc(w*h*sizeof*I0);
		for (int i = 0; i < w*h; i++) {
			I0[i] = I1[i] = 0;
			for (int l = 0; l < pd; l++) {
				I0[i] += in_a[i*pd+l];
				I1[i] += in_b[i*pd+l];
			}
		}
		void Dual_TVL1_optic_flow_multiscale(float *I0, float *I1,
			float *u1, float *u2,
			const int nxx, const int nyy,
			const float tau, const float lambda, const float theta,
			const int nscales, const float zfactor, const int warps,
			const float epsilon, const bool  verbose);
		if (global_nparams != 5)
			fail("tvl1 needs 5 pars: \"tau lbd theta nwarp eps\"");
		float p_tau = global_params[0];
		float p_lambda = global_params[1];
		float p_theta = global_params[2];
		int p_nscales = abs(global_nscales);
		float p_zfactor = 1.0/fabs(global_scale_step);
		int p_nwarps = global_params[3];
		float p_epsilon = global_params[4];
		fprintf(stderr, "tvl1 %g %g %g %d %g %d %g\n", p_tau, p_lambda,
				p_theta, p_nscales, p_zfactor, p_nwarps,
				p_epsilon);
		Dual_TVL1_optic_flow_multiscale(I0, I1, u, v, w, h,
				p_tau, p_lambda, p_theta, p_nscales, p_zfactor,
				p_nwarps, p_epsilon, 0);
		//fprintf(stderr, "tvl1 finished\n");
		free(I0);
		free(I1);
	} else {
		multi_scale_optical_flow_pd(global_method_id,
			global_params, global_nparams,
			u, v, in_a, in_b, w, h, pd,
			global_nscales, global_scale_step, global_last_scale);
	}

	switch(global_flow_visualization) {
	case VFLOW_COLORS:
		for (int i = 0; i < w*h; i++)
			flow_to_frgb(out+3*i, v[i], u[i], global_vscale);
		break;
	case VFLOW_ARROWS:
		flow_to_farr(out, u, v, w, h, global_vscale);
		break;
	case VFLOW_DIVERGENCE:
		flow_to_fdiv(out, u, v, w, h, global_vscale);
		break;
	case VFLOW_BACK:
		flow_to_fback(out, in_b, u, v, w, h);
		break;
	case VFLOW_BACKDIFF:
		flow_to_fbackdiff(out, in_a, in_b, u, v, w, h, global_vscale);
		break;
	case VFLOW_BARROWS:
		flow_to_barr(out, in_a, u, v, w, h, global_vscale);
		break;
	default: fail("unrecognized flow visualization type");
	}

	free(u);
	free(v);
}
// process one frame
static void process_tacu(float *out, float *in, int w, int h, int pd)
{
	// convert image to gray
	image_double in_gray = new_image_double(w, h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;
		float r = in[3*idx+0];
		float g = in[3*idx+1];
		float b = in[3*idx+2];
		in_gray->data[idx] = (r + g + b)/3;
	}

	// fill-in gray values (for visualization)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;
		out[3*idx+0] = in_gray->data[idx];
		out[3*idx+1] = in_gray->data[idx];
		out[3*idx+2] = in_gray->data[idx];
	}

	// computi harris-hessian
	double tic = seconds();
	ntuple_list hp = harris2(in_gray,
			global_harris_sigma,
			global_harris_k,
			global_harris_flat_th,
			global_harris_neigh
			);
	tic = seconds() - tic;
	fprintf(stderr, "harris took %g milliseconds (%g hz)\n",
			tic*1000, 1/tic);

	// plot detected keypoints
	int n[][2] = {
		{0,0},
	       	{-1,0}, {0,-1}, {0,1}, {1,0}, // 5
		{-1,-1}, {-1,+1}, {1,-1}, {1,+1} // 9
	}, nn = 9;

	for (int i = 0; i < hp->size; i++)
	{
		assert(hp->dim == 2);
		int x = hp->values[2*i+0];
		int y = hp->values[2*i+1];
		for (int p = 0; p < nn; p++)
		{
			int xx = x + n[p][0];
			int yy = y + n[p][1];
			int idx = yy*w + xx;
			if (idx < 0 || idx >= w*h) continue;
			out[3*idx + 0] = 0;
			out[3*idx + 1] = 255;
			out[3*idx + 2] = 0;
		}
		//int idx = y*w + x;
		//if (idx < 0 || idx >= w*h) continue;
		//out[3*idx + 0] = 0;
		//out[3*idx + 1] = 0;
		//out[3*idx + 2] = 255;
	}

	free_ntuple_list(hp);
	free_image_double(in_gray);

	// draw HUD
	char buf[1000];
	snprintf(buf, 1000, "sigma = %g\nk=%g\nt=%g\nn=%d",
			global_harris_sigma,
			global_harris_k,
			global_harris_flat_th,
			global_harris_neigh);
	float fg[] = {0, 255, 0};
	put_string_in_float_image(out,w,h,3, 5,5, fg, 0, &global_font, buf);
}

// process one frame
static void process_tacu_old(float *out, float *in, int w, int h, int pd)
{
	if (pd != 3) fail("bad pd");
	//for (int i = 0; i < w*h*pd; i++)
	//	out[i] = in[i];

	static int count = 0;
	static float *a = NULL;
	static float *b = NULL;
	static int ow = 0;
	static int oh = 0;

	fprintf(stderr, "pt count %d {%g}\t", count, global_vscale);
	if (global_display_img || global_display_diff)
		fprintf(stderr, "\n");

	if (count > 0 && (ow!=w || oh!=h))
		fail("do not change size, please");

	if (count == 0) { // first call
		ow = w;
		oh = h;
		a = xmalloc(3*w*h*sizeof*a);
		b = xmalloc(3*w*h*sizeof*a);
		for (int i = 0; i < w*h*pd; i++)
			a[i] = in[i];
	}
	if (count == 1) { // second call
		for (int i = 0; i < w*h*pd; i++)
			b[i] = in[i];
	}
	if (count > 1) { // rest of the calls
		float *t = b; b = a; a = t;
		for (int i = 0; i < w*h*pd; i++)
			b[i] = in[i];
		d_tacu(out, a, b, w, h, pd);
	}

	if (global_display_img && count == 100)
		iio_save_image_float_vec("/tmp/a.png", a, w, h, pd);
	if (global_display_img && count == 110)
		iio_save_image_float_vec("/tmp/b.png", a, w, h, pd);

	count += 1;
}

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

int main( int argc, char *argv[] )
{
	if (argc != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s method \"params\" step nscales lastscale vnorm\n", *argv);
	//       0  1        2        3    4       5         6
		return EXIT_FAILURE;
	}

	global_method_id = argv[1];
	char *parstring = argv[2];
	global_scale_step = atof(argv[3]);
	global_nscales = atoi(argv[4]);
	global_last_scale = atoi(argv[5]);
	global_nparams = parse_floats(global_params, 0x100, parstring);
	global_vscale = atof(argv[6]);

	global_font = *xfont8x13;
	global_font = reformat_font(global_font, UNPACKED);

	CvCapture *capture = 0;
	int accum_index = 0;
	int       key = 0;

	/* initialize camera */
	capture = cvCaptureFromCAM( atoi(argv[1]) );

	/* always check */
	if ( !capture )
		fail("could not get a capture");

	IplImage *frame = cvQueryFrame(capture);
	if (!frame) fail("did not get frame");
	int w = frame->width, W = 512;
	int h = frame->height, H = 512;
	int pd = frame->nChannels;
	int depth = frame->depth;
	if (w != 640 || h != 480 || pd != 3)
		fail("unexpected webcam size, "
				"please change some hard-coded numbers");

	fprintf(stderr, "%dx%d %d [%d]\n", w, h, pd, depth);
	//if (W > w || H > h) fail("bad crop");
	CvSize size;
	size.width = W;
	size.height = H;
	IplImage *frame_small = cvCreateImage(size, depth, pd);
	//IplImage *frame_big = cvCreateImage(size, depth, pd);
	//fprintf(stderr, "%dx%d %d [%d]\n", frame_big->width, frame_big->height, pd, depth);

	float *taccu_in = xmalloc(W*H*pd*sizeof*taccu_in);
	float *taccu_out = xmalloc(W*H*pd*sizeof*taccu_in);
	for (int i = 0; i < W*H; i++) {
		int g = 0;//rand()%0x100;
		taccu_in[3*i+0] = g;
		taccu_in[3*i+1] = g;
		taccu_in[3*i+2] = g;
	}

	/* create a window for the video */
	cvNamedWindow( "result", CV_WINDOW_FREERATIO );
	cvResizeWindow("result", W, H);

	while( key != 'q' ) {
		/* get a frame */
		frame = cvQueryFrame( capture );

		/* always check */
		if( !frame ) break;

		if (frame->width != w) fail("got bad width");
		if (frame->height != h) fail("got bad height");
		if (frame->nChannels != pd) fail("got bad nc");
		if (frame->depth != depth) fail("got bad depth");
		if (pd != 3) fail("pd is not 3");

		//for (int i = 0; i < W * H * pd; i++)
		//{
		//	taccu_in[i] = (float)(unsigned char)frame->imageData[i];
		//}
		for (int j = 0; j < 384; j++)
		for (int i = 0; i < 512; i++)
		for (int l = 0; l < pd; l++)
			taccu_in[((j+64)*512+i)*pd+l] = (float)(unsigned char)
				frame->imageData[((j+48)*w+i+64)*pd+l];

		process_tacu(taccu_out, taccu_in, W, H, pd);

		taccu_out[0]=taccu_out[1]=taccu_out[2]=0;
		taccu_out[3]=taccu_out[4]=taccu_out[5]=255;

		for (int i = 0; i < W * H * pd; i++)
			frame_small->imageData[i] = taccu_out[i];

		cvShowImage( "result", frame_small );

		/* exit if user press 'q' */
		key = cvWaitKey( 1 ) % 0x10000;
		double wheel_factor = 1.1;
		if (key == 's') global_harris_sigma /= wheel_factor;
		if (key == 'S') global_harris_sigma *= wheel_factor;
		if (key == 'k') global_harris_k /= wheel_factor;
		if (key == 'K') global_harris_k *= wheel_factor;
		if (key == 't') global_harris_flat_th /= wheel_factor;
		if (key == 'T') global_harris_flat_th *= wheel_factor;
		if (key == 'n' && global_harris_neigh > 1)
			global_harris_neigh -= 1;
		if (key == 'N') global_harris_neigh += 1;
		if (isalpha(key)) {
			printf("harris_sigma = %g\n", global_harris_sigma);
			printf("harris_k = %g\n", global_harris_k);
			printf("harris_t = %g\n", global_harris_flat_th);
			printf("harris_n = %d\n", global_harris_neigh);
			printf("\n");
		}
#if 0 // commented for surys demo
		if (key == 'd') {
			global_display_diff = !global_display_diff;
			global_display_img = false;
		}
		if (key == 'i') {
			global_display_img = !global_display_img;
			global_display_diff = false;
		}
		if (key == 'a') global_flow_visualization = VFLOW_ARROWS;
		if (key == 's') global_flow_visualization = VFLOW_BARROWS;
		if (key == 'c') global_flow_visualization = VFLOW_COLORS;
		if (key == 'y') global_flow_visualization = VFLOW_DIVERGENCE;
		if (key == 'b') global_flow_visualization = VFLOW_BACK;
		if (key == 'v') global_flow_visualization = VFLOW_BACKDIFF;
		if (key == '(') global_vscale /= global_vscale_factor;
		if (key == ')') global_vscale *= global_vscale_factor;
#endif
		//if (key > 0) {
		//	fprintf(stderr, "key = %d '%c'\n", key, key);
		//	char buf[0x100];
		//	snprintf(buf, 0x100, "/tmp/bcam_%c.png", key);
		//	iio_save_image_float_vec(buf, taccu_out, W, H, pd);
		//}
	}

	/* free memory */
	cvDestroyWindow( "result" );
	cvReleaseCapture( &capture );

	return 0;
}
