#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <omp.h>

double seconds(void)
{
	static int initialized = 0;
	static time_t first_seconds;
	struct timespec t[1];
	if (!initialized) {
		clock_gettime(CLOCK_REALTIME, t);
		first_seconds = t->tv_sec;
		initialized = 1;
	}
	clock_gettime(CLOCK_REALTIME, t);
	assert(t->tv_sec >= first_seconds);
	double r = (t->tv_sec - first_seconds) + 1e-9*t->tv_nsec;
	return r;
}

float getpixel(const float *x, int w, int h, int i, int j)
{
	if (i < 0) return 0;
	if (j < 0) return 0;
	if (i >= w) return 0;
	if (j >= h) return 0;
	return x[j*w+i];
}

void compute_laplacian_naive(float *y, float *x, int w, int h)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[w*j+i] = -4*getpixel(x, w, h, i  , j  )
			    + getpixel(x, w, h, i+1, j  )
			    + getpixel(x, w, h, i-1, j  )
			    + getpixel(x, w, h, i  , j+1)
			    + getpixel(x, w, h, i  , j-1);
}

void compute_laplacian_naive2(float *y, float *x, int w, int h)
{
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
		y[w*j+i] = -4*getpixel(x, w, h, i  , j  )
			    + getpixel(x, w, h, i+1, j  )
			    + getpixel(x, w, h, i-1, j  )
			    + getpixel(x, w, h, i  , j+1)
			    + getpixel(x, w, h, i  , j-1);
}

void compute_laplacian_parallel(float *y, const float *x, const int w, const int h)
{
#pragma omp parallel for
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		y[w*j+i] = -4*getpixel(x, w, h, i  , j  )
			    + getpixel(x, w, h, i+1, j  )
			    + getpixel(x, w, h, i-1, j  )
			    + getpixel(x, w, h, i  , j+1)
			    + getpixel(x, w, h, i  , j-1);
}

void compute_laplacian_parallel2(float *y, const float *x, const int w, const int h)
{
#pragma omp parallel for
	for (int i = 0; i < w; i++)
	for (int j = 0; j < h; j++)
		y[w*j+i] = -4*getpixel(x, w, h, i  , j  )
			    + getpixel(x, w, h, i+1, j  )
			    + getpixel(x, w, h, i-1, j  )
			    + getpixel(x, w, h, i  , j+1)
			    + getpixel(x, w, h, i  , j-1);
}

void compute_laplacian_parallel3(float *y, float *x, int w, int h)
{
	omp_set_num_threads(2);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int jfirst, jlast;

		if (tid == 0) {
			jfirst = 0;
			jlast = h/2;
		} else {
			jfirst = h/2;
			jlast = h;
		}

		for (int j = jfirst; j < jlast; j++)
			for (int i = 0; i < w; i++)
				y[w*j+i] = -4*getpixel(x, w, h, i  , j  )
					+ getpixel(x, w, h, i+1, j  )
					+ getpixel(x, w, h, i-1, j  )
					+ getpixel(x, w, h, i  , j+1)
					+ getpixel(x, w, h, i  , j-1);
	}
}

void compute_laplacian_smart(float *y, float *x, int w, int h)
{
	int i, j;
	// fill central region
	for (j = 1; j < h-1; j++)
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first row
	j = 0;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    //+ x[w*(j-1)+(i)]
			    ;

	// fill center of last row
	j = h-1;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    //+ x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first column
	i = 0;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    //+ x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of last column
	i = w-1;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    //+ x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// corner
	i = 0;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

	// corner
	i = 0;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

}

void compute_laplacian_smart2(float *y, float *x, int w, int h)
{
	int i, j;
	// fill central region
	for (i = 1; i < w-1; i++)
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first row
	j = 0;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    //+ x[w*(j-1)+(i)]
			    ;

	// fill center of last row
	j = h-1;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    //+ x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first column
	i = 0;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    //+ x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of last column
	i = w-1;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    //+ x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// corner
	i = 0;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

	// corner
	i = 0;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

}

void compute_laplacian_smart_parallel(float *y, float *x, int w, int h)
{
	int i, j;
	// fill central region
#pragma omp parallel for
	for (j = 1; j < h-1; j++)
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first row
	j = 0;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    //+ x[w*(j-1)+(i)]
			    ;

	// fill center of last row
	j = h-1;
	for (i = 1; i < w-1; i++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    //+ x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of first column
	i = 0;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    + x[w*(j)+(i+1)]
			    //+ x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// fill center of last column
	i = w-1;
	for (j = 1; j < h-1; j++)
		y[w*j+i] = -4*x[w*(j)+(i)]
			    //+ x[w*(j)+(i+1)]
			    + x[w*(j)+(i-1)]
			    + x[w*(j+1)+(i)]
			    + x[w*(j-1)+(i)]
			    ;

	// corner
	i = 0;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = 0;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		+ x[w*(j+1)+(i)]
		//+ x[w*(j-1)+(i)]
		;

	// corner
	i = w-1;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		//+ x[w*(j)+(i+1)]
		+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

	// corner
	i = 0;
	j = h-1;
	y[w*j+i] = -4*x[w*(j)+(i)]
		+ x[w*(j)+(i+1)]
		//+ x[w*(j)+(i-1)]
		//+ x[w*(j+1)+(i)]
		+ x[w*(j-1)+(i)]
		;

}

void run_several_times(float *y, float *x, int w, int h, int n,
		void (*lap)(float*,float*,int,int))
{
	for (int i = 0; i < n; i++)
		lap(y, x, w, h);
}

int main(int c, char **v)
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s w h n\n", *v);
		//                          0 1 2 3
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	int n = atoi(v[3]);

	double time_counter = seconds();

	float *x = malloc(w*h*sizeof*x);
	float *y = malloc(w*h*sizeof*y);
	float *z = malloc(w*h*sizeof*z);
	if (!x || !y || !z)
		return 2;

	for (int i = 0; i < w*h; i++)
		x[i] = 255*(rand()/(double)RAND_MAX);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME initialization = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_naive);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME naive = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_naive2);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME naive2 = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_smart);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME smart = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_smart2);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME smart2 = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_parallel);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME parallel = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_parallel2);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME parallel2 = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_parallel3);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME parallel3 = %gs\n", time_counter);

	time_counter = seconds();
	run_several_times(y, x, w, h, n, compute_laplacian_smart_parallel);
	time_counter = seconds() - time_counter;
	fprintf(stderr, "TIME smart_parallel = %gs\n", time_counter);

	free(x); free(y); free(z);
	return 0;
}
