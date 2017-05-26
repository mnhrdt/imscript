#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"

static void apply_homography(double y[2], double H[9], double x[2])
{
	double a = x[0]*H[0] + x[1]*H[1] + H[2];
	double b = x[0]*H[3] + x[1]*H[4] + H[5];
	double c = x[0]*H[6] + x[1]*H[7] + H[8];
	y[0] = a / c;
	y[1] = b / c;
}

int main_nodelta(int c, char *v[])
{
	if (c != 5 && c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s w h delta [ofield] < crophoms\n", *v);
		//                0 1 2 3      5
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	//double delta = atof(v[3]);
	char *filename_out = c > 4 ? v[4] : "-";

	float *flow = malloc(2*w*h*sizeof*flow);
	int *count = malloc(w*h*sizeof*count);
	for (int i = 0; i < w*h; i++)
		flow[i] = flow[w*h+i] = count[i] = 0;

	int k[4];
	double H[9];
	int cx = 0;
	while (4+9 == scanf("%d %d %d %d\t"
				"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			k, k+1, k+2, k+3,
			H, H+1, H+2, H+3, H+4, H+5, H+6, H+7, H+8
			))
	{
		fprintf(stderr, "CX=%d k=%d %d %d %d ht=%g %g\n",
				cx++, k[0], k[1], k[2], k[3],
				H[2], H[5]);
		for (int j = k[1]; j <= k[3]; j++)
		for (int i = k[0]; i <= k[2]; i++)
			if (i >= 0 && j >= 0 && i < w && j < h)
			{
				int idx = w*j + i;
				double p[2] = {i, j}, q[2];
				apply_homography(q, H, p);
				flow[2*idx+0] += q[0] - p[0];
				flow[2*idx+1] += q[1] - p[1];
				count[idx] += 1;
			}
	}
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = w*j + i;
		if (count[idx] > 0)
		{
			flow[2*idx+0] /= count[idx];
			flow[2*idx+1] /= count[idx];
		}
	}

	iio_save_image_float_vec(filename_out, flow, w, h, 2);
	iio_save_image_int("/tmp/counts.tiff", count, w, h);

	free(count);
	free(flow);
	return 0;
}

static double dist2(double a[2], double b[2])
{
	double x = a[0] - b[0];
	double y = a[1] - b[1];
	return x*x + y*y;
}

int main(int c, char *v[])
{
	if (c != 5 && c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s w h delta [ofield] < crophoms\n", *v);
		//                0 1 2 3      5
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	double delta = atof(v[3]);
	char *filename_out = c > 4 ? v[4] : "-";

	float *flow = malloc(2*w*h*sizeof*flow);
	float *weig = malloc(w*h*sizeof*weig);
	//int *count = malloc(w*h*sizeof*count);
	for (int i = 0; i < w*h; i++)
		flow[i] = flow[w*h+i] = weig[i] = 0;

	int k[4];
	double H[9];
	int cx = 0;
	while (4+9 == scanf("%d %d %d %d\t"
				"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			k, k+1, k+2, k+3,
			H, H+1, H+2, H+3, H+4, H+5, H+6, H+7, H+8
			))
	{
		fprintf(stderr, "CX=%d k=%d %d %d %d ht=%g %g\n",
				cx++, k[0], k[1], k[2], k[3],
				H[2], H[5]);
		double center[2] = {(k[0]+k[2])/2, (k[1]+k[3])/2};
		//for (int j = k[1]; j <= k[3]; j++)
		//for (int i = k[0]; i <= k[2]; i++)
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
			//if (i >= 0 && j >= 0 && i < w && j < h)
			{
				int idx = w*j + i;
				double p[2] = {i, j}, q[2];
				apply_homography(q, H, p);
				double omega = dist2(center,p)/(2*delta*delta);
				omega = exp(-omega);
				flow[2*idx+0] += (q[0] - p[0])*omega;
				flow[2*idx+1] += (q[1] - p[1])*omega;
				weig[idx] += omega;
			}
	}
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = w*j + i;
		if (weig[idx] > 0)
		{
			flow[2*idx+0] /= weig[idx];
			flow[2*idx+1] /= weig[idx];
		}
	}

	iio_save_image_float_vec(filename_out, flow, w, h, 2);
	iio_save_image_float("/tmp/weig.tiff", weig, w, h);

	free(weig);
	free(flow);
	return 0;
}
