#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double ell_pee_distance(float *x, float *y, int n, float p)
{
	double r = 0;
	if (isinf(p)) {
		r = -INFINITY;
		for (int i = 0; i < n; i++) {
			double t = fabs(x[i] - y[i]);
			if (t > r)
				r = t;
		}
	} else if (fabs(p) == 0) {
		for (int i = 0; i < n; i++)
			if (x[i] != y[i])
				r += 1;
		if (signbit(p))
			r /= n;
	} else if (isnormal(p)) {
		for (int i = 0; i < n; i++)
			r += pow(fabs(x[i] - y[i]), fabs(p));
		r = pow(r, 1/fabs(p));
		if (p < 0)
			r /= pow(n, 1/fabs(p));
	}
	return r;

}

static double mean_square_error(float *x, float *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += (x[i] - y[i]) * (x[i] - y[i]);
	r /= n;
	return r;
}

static double root_mean_square_error(float *x, float *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r = hypot(r, x[i] - y[i]);
	r /= sqrt(n);
	return r;
}

static double mean_absolute_error(float *x, float *y, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += fabs(x[i] - y[i]);
	r /= n;
	return r;
}

static void meanvar(double *mean, double *var, float *x, int n)
{
	*mean = 0;
	for (int i = 0; i < n; i++)
		*mean += x[i];
	*mean /= n;
	for (int i = 0; i < n; i++)
		*var += (x[i] - *mean) * (x[i] - *mean);
	*var /= n - 1;
}

static void minmax(double *min, double *max, float *x, int n)
{
	if (min) {
		*min = INFINITY;
		for (int i = 0; i < n; i++)
			if (x[i] < *min)
				*min = x[i];
	}
	if (max) {
		*max = -INFINITY;
		for (int i = 0; i < n; i++)
			if (x[i] > *max)
				*max = x[i];
	}
}

static double dynamic_range(float *x, int n)
{
	double min, max;
	minmax(&min, &max, x, n);
	return max - min;
}


static double uiqi(float *x, float *y, int n)
{
	double mx, my, sx, sy, sxy = 0;
	meanvar(&mx, &sx, x, n);
	meanvar(&my, &sy, y, n);
	for (int i = 0; i < n; i++)
		sxy += (x[i] - mx)*(y[i] - my);
	sxy /= n - 1;
	double q = (4 * sxy * mx * my);
	q /= (sx + sy) * (mx*mx + my*my);
	return q;
}


#define SSIM_K1 0.01
#define SSIM_K2 0.03

static double ssim(float *x, float *y, int n)
{
	double mx, my, sx, sy, sxy = 0;
	meanvar(&mx, &sx, x, n);
	meanvar(&my, &sy, y, n);
	for (int i = 0; i < n; i++)
		sxy += (x[i] - mx)*(y[i] - my);
	sxy /= n - 1;
	double Lx = dynamic_range(x, n);
	double Ly = dynamic_range(x, n);
	double L = (Lx + Ly)/2;
	double C1 = SSIM_K1 * L;
	double C2 = SSIM_K2 * L;
	C1 *= C1;
	C2 *= C2;
	double r = (2*mx*my + C1) * (2*sxy + C2);
	r /= (mx*mx + my*my + C1) * (sx + sx + C2);
	return r;
}

static double psnr(float *x, float *y, int n)
{
	return 20 * log10(dynamic_range(x, n)/root_mean_square_error(x, y, n));
}

static bool string_is_lp(char *s, double *p)
{
	if (s[0] == 'L') {
		double pee;
		int r = sscanf(s+1, "%lf", &pee);
		if (r == 1 && !isnan(pee)) {
			*p = pee;
			return true;
		}
	}
	return false;
}

static double ncc(float *x, float *y, int n)
{
	double mx = 0, my = 0;
	for (int i = 0; i < n; i++) {
		mx += x[i];
		my += y[i];
	}
	mx /= n;
	my /= n;

	double sx = 0, sy = 0, sxy = 0;
	for (int i = 0; i < n; i++) {
		sx  += (x[i] - mx) * (x[i] - mx);
		sy  += (y[i] - my) * (y[i] - my);
		sxy += (x[i] - mx) * (y[i] - my);
	}
	sx = sqrt(sx/n);
	sy = sqrt(sy/n);

	return sxy/(n*sx*sy);
}

#include "fail.c"

static double imgerr(char *m, float *x, float *y, int n)
{
	double r;
	if (false);
	else if (0 == strcmp(m, "MSE"))   r = mean_square_error(x, y, n);
	else if (0 == strcmp(m, "RMSE"))  r = root_mean_square_error(x, y, n);
	else if (0 == strcmp(m, "MAE"))   r = mean_absolute_error(x, y, n);
	else if (0 == strcmp(m, "UIQI"))  r = uiqi(x, y, n);
	else if (0 == strcmp(m, "SSIM"))  r = ssim(x, y, n);
	else if (0 == strcmp(m, "PSNR"))  r = psnr(x, y, n);
	else if (0 == strcmp(m, "NCC"))   r = ncc(x, y, n);
	else if (string_is_lp(m, &r))     r = ell_pee_distance(x, y, n, r);
	else fail("unrecognized metric \"%s\"", m);
	return r;
}


#include "iio.h"

#include "smapa.h"
#include "xmalloc.c"
SMART_PARAMETER_SILENT(IMGERR_PPT,0)

int main(int c, char *v[])
{
	if (c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s metric imga [imgb]\n", *v);
		//                          0 1      2     3
		return 1;
	}
	char *metric_id = v[1];
	char *filename_a = v[2];
	char *filename_b = c > 3 ? v[3] : "-";

	int w[2], h[2], pd[2];
	float *x = iio_read_image_float_vec(filename_a, w+0, h+0, pd+0);
	float *y = iio_read_image_float_vec(filename_b, w+1, h+1, pd+1);
	if (w[0] != w[1] || h[0] != h[1] || pd[0] != pd[1])
		fail("input image sizes mismatch (%g %g %g) != (%g %g %g)",
					w[0], h[0], pd[0], w[1], h[1], pd[1]);
	float *xx = x, *yy = y;
	int m = IMGERR_PPT();
	if (m > 0 && 2*m+1<*w && 2*m+1<*h) {
		xx = xmalloc((*w-2*m)*(*h-2*m)**pd*sizeof*xx);
		yy = xmalloc((*w-2*m)*(*h-2*m)**pd*sizeof*yy);
		for (int i = 0; i < *w-2*m; i++)
		for (int j = 0; j < *h-2*m; j++)
		for (int l = 0; l < *pd; l++)
		{
			xx[((*w-2*m)*j + i)**pd+l] = x[(*w*(j+m)+(i+m))**pd+l];
			yy[((*w-2*m)*j + i)**pd+l] = y[(*w*(j+m)+(i+m))**pd+l];
		}
		*w -= 2*m;
		*h -= 2*m;
	}

	iio_write_image_float_vec("/tmp/xx.tmp.asc", xx, *w, *h, *pd);
	iio_write_image_float_vec("/tmp/yy.tmp.asc", yy, *w, *h, *pd);
	double e = imgerr(metric_id, xx, yy, *w**h**pd);

	printf("%.16lf\n", e);
	return 0;
}
