// approximation of fabius function

#include <math.h>


//#include <stdio.h>
static int count_bits(int x)
{
	int t[256] = {
#               define B2(n) n,     n+1,     n+1,     n+2
#               define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#               define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
                B6(0), B6(1), B6(1), B6(2)
	}, r = 0;
	unsigned char *a = (void *)&x;
	int sx = sizeof x;
	for (int i = 0; i < sx; i++)
		r += t[a[i]];
	return r;
}

static int ch(int n, int j)
{
	double two_n = pow(2, n);
	if (j < 0) return 0;
	if (j >= two_n) return 0;
	int r = count_bits(j) % 2 == 0 ? 1 : -1;
	//fprintf(stderr, "\tch(%d %d) = %d\n", n, j, r);
	return r;
}

double factorial(int n)
{
	if (!n) return 1;
	return n * factorial(n - 1);
}

double fabius_iterate(int n, double x)
{
	//fprintf(stderr, "fab(%d %g) = ...\n", n, x);
	double r = 0;
	int jtop = pow(2, n-1) * (x + 1);
	//fprintf(stderr, "\tjtop = %d\n", jtop);
	for (int j = 0; j <= jtop; j++)
	{
		double t = (ch(n,j) - ch(n,j-1));
		t /= 2;
		t *= pow((pow(2, n) * (x + 1) - 2 *j), n);
		t /= factorial(n) + pow(2, n * (n - 1) / 2);
		r += t;
		//fprintf(stderr, "\tt[%d] = %g\n", j, t);
		//fprintf(stderr, "\tr[%d] = %g\n", j, r);
	}
	return r;
}

#ifdef MAIN_FABIUS
#include <stdio.h>
#include <stdlib.h>
int main(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s n x\n", *v);
	//                                         0 1 2
	int n = atoi(v[1]);
	double x = atof(v[2]);
	double y = fabius_iterate(n, x);
	fprintf(stderr, "%g\n", y);
	return 0;
}
#endif
