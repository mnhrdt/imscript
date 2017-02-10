
#ifndef ODDP
#define ODDP(x) ((x)&1)
#endif
#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

#include "fail.c"
#include "xmalloc.c"

struct statistics_float {
	float min, max, median, average, sample, variance, middle;
};

void print_stats(FILE *f, struct statistics_float *s, char *name)
{
#ifdef _STDIO_H
	fprintf(f, "%s min: %g\n", name, s->min);
	fprintf(f, "%s max: %g\n", name, s->max);
	fprintf(f, "%s median: %g\n", name, s->median);
	fprintf(f, "%s average: %g\n", name, s->average);
	fprintf(f, "%s variance: %g\n", name, s->variance);
#endif//STDIO_H
}

//SMART_PARAMETER_INT(STATISTIC_MEDIAN_BIAS,0)
//SMART_PARAMETER_SILENT(STATISTIC_MIDDLE_BIAS,-1)

#define STATISTIC_MEDIAN_BIAS 0
#define STATISTIC_MIDDLE_BIAS -1

static int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void statistics_getf_spoilable(struct statistics_float *s, float *f,
		int n)
{
	s->middle = f[n/2-1];
	int mt = STATISTIC_MEDIAN_BIAS;
	int mi = STATISTIC_MIDDLE_BIAS;
	switch(mi)
	{
		case -1: break;
		case 0: s->middle += f[n/2]; s->middle /=2; break;
		case 1: s->middle = f[n/2]; break;
		default: fail("bad STATISTIC_MEDIAN_BIAS %d", mt);
	}
	//
	qsort(f, n, sizeof*f, compare_floats);
	s->min = f[0];
	s->max = f[n-1];
	s->median = f[n/2-1];
	if (EVENP(n))
	{
		int mtype = STATISTIC_MEDIAN_BIAS;
		switch(mtype)
		{
			case -1: break;
			case 0: s->median += f[n/2]; s->median /=2; break;
			case 1: s->median = f[n/2]; break;
			default: fail("bad STATISTIC_MEDIAN_BIAS %d", mtype);
		}
	}
	s->average = 0;
	for (int i = 0; i < n; i++)
		s->average += f[i];
	s->average /= n;
	s->variance = 0;
	for (int i = 0; i < n; i++)
		s->variance = hypot(s->variance, s->average - f[i]);
	s->sample = f[randombounds(0, n-1)];
}

static void statistics_getf(struct statistics_float *s, float *fin, int n)
{
	if (n > 1000) {
		float *f = xmalloc(n*sizeof*f);
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
		xfree(f);
	} else {
		float f[n];
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
	}
}

void statistics_getf_stride(struct statistics_float *s,
				float *fin, int n, int stride)
{
	float *f = xmalloc(n*sizeof*f);
	for(int i = 0; i < n; i++)
		f[i] = fin[stride*i];
	statistics_getf_spoilable(s, fin, n);
	xfree(f);
}

