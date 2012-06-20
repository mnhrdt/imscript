#ifndef _RANDOM_C
#define _RANDOM_C

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


static double random_uniform(void)
{
	return rand()/(RAND_MAX+1.0);
}

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y1 = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y1;
}

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

#endif//_RANDOM_C
