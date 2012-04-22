/**
 * Display video from webcam
 *
 * Author  Nash
 * License GPL
 * Website http://nashruddin.com
 */

#include <stdio.h>
#include "cv.h"
#include "highgui.h"

#define MAX_FRAME_ACCUM 100

static void stats(float *x, int n, int d)
{
	float mi = x[0];
	float ma = x[0];
	for (int i = 1; i < n*d; i++)
	{
		if (x[i] < mi) mi = x[i];
		if (x[i] > ma) ma = x[i];
	}
	fprintf(stderr, "[%g %g]\n", mi, ma);
}


#include "xmalloc.c"

int main( int argc, char **argv )
{
	if (argc != 3) return 1;
	int accum_n = atoi(argv[2]);
	if (accum_n > MAX_FRAME_ACCUM)
		accum_n = MAX_FRAME_ACCUM;
	if (accum_n < 1)
		accum_n = 1;


	CvCapture *capture = 0;
	float *accum[MAX_FRAME_ACCUM+1];
	int accum_index = 0;
	int       key = 0;
	int count = 0;

	/* initialize camera */
	capture = cvCaptureFromCAM( atoi(argv[1]) );

	/* always check */
	if ( !capture )
		fail("could not get a capture");

	IplImage *frame = cvQueryFrame(capture);
	if (!frame) fail("did not get frame");
	int w = frame->width;
	int h = frame->height;
	int pd = frame->nChannels;
	int depth = frame->depth;
	for (int i = 0; i < accum_n + 1; i++) {
		accum[i] = xmalloc(w * h * pd * sizeof(float));
		for (int j = 0; j < w*h*pd; j++)
			accum[i][j] = 0;
	}

	fprintf(stderr, "%dx%d %d [%d]\n", w, h, pd, depth);
	CvSize size;
	size.width = 2*w;
	size.height = 2*h;
	IplImage *frame_big = cvCreateImage(size, depth, pd);
	fprintf(stderr, "%dx%d %d [%d]\n", frame_big->width, frame_big->height, pd, depth);

	/* create a window for the video */
	cvNamedWindow( "result", CV_WINDOW_FREERATIO );
	//cvResizeWindow("result", 2*w, 2*h);
	cvResizeWindow("result", w, h);

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

		// copy this frame values into the carroussel
		for (int i = 0; i < w * h * pd; i++)
			accum[accum_index][i] = (unsigned char)frame->imageData[i];

		//fprintf(stderr, "accum_index = %d\n", accum_index);

		float *taccu = accum[accum_n];
		for (int i = 0; i < w * h * pd; i++)
			taccu[i] = 0;//(float)(unsigned char)frame->imageData[i];


		for (int j = 0; j < accum_n; j++)
			for (int i = 0; i < w * h * pd; i++)
				taccu[i] += accum[j][i]/(accum_n);
		//
		//if (0)if (count++ == 9 && pd == 3)
		//	for (int i = 0; i < w * h * 3; i++)
		//		fprintf(stderr, "%d %d %d\n",
		//				frame->imageData[3*i],
		//				frame->imageData[3*i+1],
		//				frame->imageData[3*i+2]);

		//stats(taccu,w*h,pd);
		taccu[0]=taccu[1]=taccu[2]=0;
		taccu[3]=taccu[4]=taccu[5]=255;

		for (int i = 0; i < (w * h * pd)/2; i++)
			frame->imageData[i] = taccu[i];

		//for (int i = 0; i < w * h * pd; i++)
		//	frame->imageData[i] *= 4;//taccu[i];

		//for (int j = 0; j < h; j++)
		//for (int i = 0; i < w; i++)
		//for (int l = 0; l < pd; l++)
		//{
		//	int g = frame->imageData[(i+j*w)*pd+l];
		//	frame_big->imageData[((2*i+0)+(2*j+0)*2*w)*pd+l] = g;
		//	frame_big->imageData[((2*i+0)+(2*j+1)*2*w)*pd+l] = g;
		//	frame_big->imageData[((2*i+1)+(2*j+0)*2*w)*pd+l] = g;
		//	frame_big->imageData[((2*i+1)+(2*j+1)*2*w)*pd+l] = g;
		//}

		/* display current frame */
		//cvShowImage( "result", frame_big );
		cvShowImage( "result", frame );

		/* exit if user press 'q' */
		key = cvWaitKey( 1 );
		accum_index = (accum_index + 1)%accum_n;
	}

	/* free memory */
	cvDestroyWindow( "result" );
	cvReleaseCapture( &capture );

	return 0;
}
