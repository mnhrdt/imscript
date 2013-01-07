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
	{
		double t = x[i] - y[i];
		r += t*t;
	}
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
	{
		double t = x[i] - *mean;
		*var += t*t;
	}
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
	else if (string_is_lp(m, &r))     r = ell_pee_distance(x, y, n, r);
	else fail("unrecognized metric \"%s\"", m);
	return r;
}


#include "iio.h"

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

	double e = imgerr(metric_id, x, y, *w**h**pd);

	printf("%g\n", e);
	return 0;
}
