/**
 * Display video from webcam
 *
 * Author  Nash
 * License GPL
 * Website http://nashruddin.com
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "cv.h"
#include "highgui.h"



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

void multi_scale_optical_flow_pd(char *algorithm_name, float *pars, int npars,
		float *u, float *v, float *a, float *b, int w, int h, int pd,
		int nscales, float scale_step, int last_scale);

// visualize the "difference" between two consecutive frames
static void d_tacu(float *out, float *in_a, float *in_b, int w, int h, int pd)
{
	//for (int i = 0; i < w*h*pd; i++) {
	//	float g = (in_b[i] - in_a[i])/2 + 128;
	//	if (g > 255) g = 255;
	//	if (g < 0) g = 0;
	//	out[i] = g;
	//}
	//return;

	float *u = xmalloc(w*h*sizeof*u);
	float *v = xmalloc(w*h*sizeof*u);

	float params[3] = {19, 5, 0.01};
	multi_scale_optical_flow_pd("hs", params, 3,
			u, v, in_a, in_b, w, h, pd, 3, -2, 3);

	for (int i = 0; i < w*h; i++)
	{
		flow_to_frgb(out+3*i, u[i], v[i], 1);
	}

	free(u);
	free(v);
}

// process one frame
static void process_tacu(float *out, float *in, int w, int h, int pd)
{
	if (pd != 3) fail("bad pd");
	//for (int i = 0; i < w*h*pd; i++)
	//	out[i] = in[i];

	static int count = 0;
	static float *a = NULL;
	static float *b = NULL;
	static int ow = 0;
	static int oh = 0;

	fprintf(stderr, "pt count %d\n", count);

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

	count += 1;
}

#include "iio.h"
int main( int argc, char **argv )
{

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
		int g = rand()%0x100;
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
		key = cvWaitKey( 1 );
		if (key > 0) {
			fprintf(stderr, "key = %d '%c'\n", key, key);
			char buf[0x100];
			snprintf(buf, 0x100, "/tmp/bcam_%c.png", key);
			iio_save_image_float_vec(buf, taccu_out, W, H, pd);
		}
	}

	/* free memory */
	cvDestroyWindow( "result" );
	cvReleaseCapture( &capture );

	return 0;
}
