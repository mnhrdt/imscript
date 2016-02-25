#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "fail.c"
#include "bicubic.c"
#include "bilinear_interpolation.c"

static float ipos(int i, int n)
{
	return 2*3.1416*(i - n/2.0)/n;
}

static void cline(float *musigma,
		float *l, int n, float *x, int w, int h, float angle)
{
	float an = angle*3.1416/180;
	float c = cos(an);
	float s = sin(an);
	float t[2] = {(fmin(w,h)-1)/2.0, 0};
	float dir[2] = {c*t[0]+s*t[1], -s*t[0]+c*t[1]};
	float zer[2] = {(w-1.0)/2, (h-1.0)/2};
	float from[2], toto[2];
	for (int i = 0; i < 2; i++) {
		from[i] = zer[i] - dir[i];
		toto[i] = zer[i] + dir[i];
	}
	fprintf(stderr, "c,s = %g %g\n", c, s);
	fprintf(stderr, "dir = %g %g\n", dir[0], dir[1]);
	fprintf(stderr, "zer = %g %g\n", zer[0], zer[1]);
	fprintf(stderr, "(%g %g)=>(%g %g)\n",from[0],from[1],toto[0],toto[1]);
	double mu = 0, sigma = 0, nn = 0;
	for (int i = 0; i < n; i++) {
		float a = i/(n - 1.0), p[2];
		for (int j = 0; j < 2; j++)
			p[j] = (1-a)*from[j] + a*toto[j];
		//bicubic_interpolation(l + i, x, w, h, 1, p[0], p[1]);
		bilinear_interpolation_vec_at(l + i, x, w, h, 1, p[0], p[1]);
		if (isfinite(l[i])) {
			//double posi = hypot(p[0]-zer[0], p[1]-zer[1]);
			//if (p[0] < 0) posi = -posi;
			double posi = ipos(i, n);
			mu += l[i]*posi;
			nn += l[i];
		}
	}
	mu /= nn;
	for (int i = 0; i < n; i++) {
		float a = i/(n - 1.0), p[2];
		for (int j = 0; j < 2; j++)
			p[j] = (1-a)*from[j] + a*toto[j];
		double posi = ipos(i, n);
		if (isfinite(l[i]))
			sigma += l[i]*(posi-mu)*(posi-mu);
	}
	sigma = sqrt(sigma/nn);
	if (musigma) {
		musigma[0] = mu;
		musigma[1] = sigma;
	}
	fprintf(stderr, "mu sigma = %g %g\n", mu, sigma);
}

static void clineh(float *musigma, float *l, int n, float *x, int w, int h)
{
	//if (w != h) fail("clineh only works over square images");
	double mu = 0, sigma = 0, nn = 0;
	int lc[n];
	for (int i = 0; i < n; i++)
		lc[i] = l[i] = 0;
	float zer[2] = {(w-1.0)/2, (h-1.0)/2};
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float rp[2] = {(i - zer[0])/zer[0], (j - zer[1])/zer[1]};
		float rpn = hypot(rp[0], rp[1])/sqrt(2);
		float rpa = 180*atan2(rp[1], rp[0])/3.1416;
		if (rpa < 0) rpn = -rpn;
		//fprintf(stderr, "(%d %d) => (%g %g) => (%g _ %g)\n",
		//		i, j, rp[0], rp[1], rpn, rpa);
		int lidx = n*(rpn + 1)/2;
		if (lidx < 0) lidx = 0;
		if (lidx >= n) lidx = n-1;
		lc[lidx] += 1;
		l[lidx] += x[w*j+i];
	}
	for (int i = 0; i < n; i++)
		if (lc[i] > 0)
			l[i] /= lc[i];
	if (musigma) {
		musigma[0] = mu;
		musigma[1] = sigma;
	}
}

static void getmusigmamass(float musigmamass[3], float *l, int n)
{
	double mu = 0, sigma = 0, mass = 0, nn = 0;
	for (int i = 0; i < n; i++)
		if (isfinite(l[i]))
			mass += l[i];
	mass *= 2*3.1416/n;
	for (int i = 0; i < n; i++)
		if (isfinite(l[i])) {
			mu += l[i] * ipos(i, n);
			nn += l[i];
		}
	mu /= nn;
	for (int i = 0; i < n; i++) {
		double posi = ipos(i, n);
		if (isfinite(l[i]))
			sigma += l[i]*(posi-mu)*(posi-mu);
	}
	sigma = sqrt(sigma/nn);
	fprintf(stderr, "mu sigma mass nn = %g %g %g %g\n", mu, sigma, mass,nn);
	musigmamass[0] = mu;
	musigmamass[1] = sigma;
	musigmamass[2] = mass;
}

static void plot_cline2(float *l, int n)
{
	printf("set samples 1000\n");
	//printf("set logscale x\n");
	//printf("set logscale y\n");
	printf("plot \"-\" w lines title \"data\"");
	float musigmamass[3];
	getmusigmamass(musigmamass, l, n);
	float mu = musigmamass[0];
	float sigma = musigmamass[1];
	float mass = musigmamass[2];
	float alpha_g = mass/(sigma*sqrt(2*3.1416));
	float alpha_l = mass/(sigma*sqrt(2));
	fprintf(stderr, "mass = %g\n", mass);
	//printf(",(%g)*exp(-(x-(%g))**2/(2*(%g)**2)) title \"normal\"",
	//		alpha_g, mu, sigma);
	//printf(",(%g)*exp(-abs(x-(%g))*sqrt(2)/(%g)) title \"laplacian\"",
	//		alpha_l, mu, sigma);
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		float ipoz = ipos(i, n);
		printf("\t%g %g\n", ipoz, l[i]);
	}
	printf("end\n");
}

static void plot_cline(float *l, int n, char *title, float mu, float sigma)
{
	if (title)
		printf("set title \"%s\"\n", title);
	printf("set samples 1000\n");
	printf("plot \"-\" w lines title \"data\"");
	if (isfinite(mu) && sigma > 0) {
		double mass = 0;
		for (int i = 0; i < n; i++)
			if (isfinite(l[i]))
				mass += l[i];
		mass *= 2*3.1416/n;
		float alpha_g = mass/(sigma*sqrt(2*3.1416));
		float alpha_l = mass/(sigma*sqrt(2));
		fprintf(stderr, "mass = %g\n", mass);
		printf(",%g*exp(-(x-%g)**2/(2*(%g)**2)) title \"normal\"",
				alpha_g, mu, sigma);
		printf(",%g*exp(-abs(x-%g)*sqrt(2)/%g) title \"laplacian\"",
				alpha_l, mu, sigma);
	}
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		float ipoz = ipos(i, n);
		printf("\t%g %g\n", ipoz, l[i]);
	}
	printf("end\n");
}

#include "iio.h"
#include "smapa.h"
SMART_PARAMETER(NFAC,1)
int main(int c, char *v[])
{
	if (c != 2 && c != 3) {
		fprintf(stderr, "usage:\n\t%s angle [img] >plot\n", *v);
		//                          0 1      2
		return 1;
	}
	float angle = atof(v[1]);
	char *in_img = c>2 ? v[2] : "-";
	int w, h;
	float *x = iio_read_image_float(in_img, &w, &h);
	int n = NFAC()*(2+hypot(w+2,h+2));
	fprintf(stderr, "geometry = %dx%d\n", w, h);
	fprintf(stderr, "n = %d\n", n);
	float l[n], musigma[2] = {0};
	if (isfinite(angle))
		cline(musigma, l, n, x, w, h, angle);
	else
		clineh(NULL, l, n, x, w, h);
	//plot_cline(l, n, v[1], musigma[0], musigma[1]);
	plot_cline2(l, n);
	free(x);
	return 0;
}
