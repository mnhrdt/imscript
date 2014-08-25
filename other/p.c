//typedef void (*f_t)(float *);
//
//void f(float *x) {}
//void g(int *x) {}
//
//void m(void)
//{
//	f_t a, b;
//	a = f;
//	b = g;
//}

typedef int (*f_t)(double);

extern f_t *f;

int f(double x)
{
	return (int)x;
}

int g(float x)
{
	return (int)x;
}

int h(float x, int y)
{
	return y+(int)x;
}

//f_t ff = f;
//f_t gg = g;
//f_t hh = h;
