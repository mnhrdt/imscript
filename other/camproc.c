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

#include "fail.c"
#include "xmalloc.c"
#include "getpixel.c"

static float getsample_avg(float *x, int w, int h, int pd, int i, int j, int L)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	float acc = 0;
	for (int l = 0; l < pd; l++)
	{
		acc += getsample_0(x, w, h, pd, i, j, l);
		acc += getsample_0(x, w, h, pd, i+1, j, l);
		acc += getsample_0(x, w, h, pd, i, j+1, l);
		acc += getsample_0(x, w, h, pd, i+1, j+1, l);
	}
	return acc/(4*pd);
}
static void inplace_blur(float *in, int w, int h, int pd)
{
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	for (int l = 0; l < pd; l++)
	{
		int idxa = (j*w+i)*pd + l;
		int idxb = (j*w+i+1)*pd + l;
		int idxc = (j*w+i+w)*pd + l;
		int idxd = (j*w+i+w+1)*pd + l;
		in[idxa] = (in[idxa] + in[idxb] + in[idxc] + in[idxd])/4;
	}
}

// process one frame
static void process_tacu(float *out, float *in, int w, int h, int pd)
{
	//for (int i = 0; i < w*h*pd; i++) out[i] = in[i];
	//return;
	getsample_operator p = getsample_0;

	for (int k = 0; k < 20; k++)
		inplace_blur(in, w, h, pd);

	if (pd != 3) fail("bad pd");
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
	{
		//int l = 0;
		int idx = (j*w+i)*pd + l;
		float p00 = p(in,w,h,pd, i, j, l);
		float p10 = p(in,w,h,pd, i+1, j, l);
		float pm0 = p(in,w,h,pd, i-1, j, l);
		float p01 = p(in,w,h,pd, i, j+1, l);
		float p0m = p(in,w,h,pd, i, j-1, l);
		float p11 = p(in,w,h,pd, i+1, j+1, l);
		float pm1 = p(in,w,h,pd, i-1, j+1, l);
		float pmm = p(in,w,h,pd, i-1, j-1, l);
		float p1m = p(in,w,h,pd, i+1, j-1, l);
		float dx = p10 - p00;
		float dy = p01 - p00;
		float lap = 4 * p00 - p10 - pm0 - p01 - p0m;
		float dxx = pm0 + p10 - 2 * p00;
		float dyy = p0m + p01 - 2 * p00;
		//float dxy = (p11 + pmm - p1m - pm1)/4;
		float dxy = (2*p00+p11+pmm-p10-p01-pm0-p0m)/2;
		float trh = dxx + dyy;//same as lap
		float doh = dxx*dyy - dxy*dxy;
		//float lv = 4*hypot(dx,dy);
		float lv = 127+17*lap;
		//float lv = 255-(0*127+4*fmax(0,doh));
		if (lv < 0) lv = 0;
		if (lv >255) lv = 255;
		out[idx+l] = lv;//(dxx * dyy - dxy * dxy)/94;
		//out[idx+1] = lv;//(dxx * dyy - dxy * dxy)/94;
		//out[idx+2] = lv;//(dxx * dyy - dxy * dxy)/94;
		//out[idx] = 127 + (dxx * dyy - dxy * dxy)/94;
		//out[idx] = p00;//127 + (dxx * dyy - dxy * dxy)/94;
	}
	return;

//	static int count = 0;
//	static float *a = NULL;
//	static float *b = NULL;
//	static int ow = 0;
//	static int oh = 0;
//
//	fprintf(stderr, "pt count %d {%g}\t", count, global_vscale);
//	if (global_display_img || global_display_diff)
//		fprintf(stderr, "\n");
//
//	if (count > 0 && (ow!=w || oh!=h))
//		fail("do not change size, please");
//
//	if (count == 0) { // first call
//		ow = w;
//		oh = h;
//		a = xmalloc(3*w*h*sizeof*a);
//		b = xmalloc(3*w*h*sizeof*a);
//		for (int i = 0; i < w*h*pd; i++)
//			a[i] = in[i];
//	}
//	if (count == 1) { // second call
//		for (int i = 0; i < w*h*pd; i++)
//			b[i] = in[i];
//	}
//	if (count > 1) { // rest of the calls
//		float *t = b; b = a; a = t;
//		for (int i = 0; i < w*h*pd; i++)
//			b[i] = in[i];
//		d_tacu(out, a, b, w, h, pd);
//	}
//
//	if (global_display_img && count == 100)
//		iio_save_image_float_vec("/tmp/a.png", a, w, h, pd);
//	if (global_display_img && count == 110)
//		iio_save_image_float_vec("/tmp/b.png", a, w, h, pd);
//
//	count += 1;
}

int main( int argc, char *argv[] )
{
	if (argc != 1)
		return 1;

	CvCapture *capture = 0;

	/* initialize camera */
	capture = cvCaptureFromCAM( 0 );
	int key = 0;

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
		key = cvWaitKey( 1 ) % 0x10000;
		//if (key == 'd') {
		//	global_display_diff = !global_display_diff;
		//	global_display_img = false;
		//}
		//if (key == 'i') {
		//	global_display_img = !global_display_img;
		//	global_display_diff = false;
		//}
		//if (key == 'a') global_flow_visualization = VFLOW_ARROWS;
		//if (key == 's') global_flow_visualization = VFLOW_BARROWS;
		//if (key == 'c') global_flow_visualization = VFLOW_COLORS;
		//if (key == 'y') global_flow_visualization = VFLOW_DIVERGENCE;
		//if (key == 'b') global_flow_visualization = VFLOW_BACK;
		//if (key == 'v') global_flow_visualization = VFLOW_BACKDIFF;
		//if (key == '(')
		//	global_vscale /= global_vscale_factor;
		//if (key == ')')
		//	global_vscale *= global_vscale_factor;
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
