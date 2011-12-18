#include <assert.h>
#include <math.h> // for fabs

static
float march_regular(float a, float b, float c, float d, float x, float y)
{
	assert(a <= b); assert(b <= c); assert(c <= d);

	if (x*(c - a) <= (1-y)*(c - a)) // first triangle
		return (c - a)*x + (b - a)*y + a;
	if ((d - b)*(x - 1) <= (d - c)*y) // second triangle
		return (d - b)*(x - 1) + (d - c)*y + c;

	// TODO: simplify formula below
	float alpha = c == a ? 0.5 : (b - a)/(c - a);
	float beta = d == b ? 0.5 : (c - b)/(d - b);
	float x0 = alpha*(1 - y);
	float x1 = 1 - (1-beta)*y;
	assert(x0 >= 0); assert(x0 <= 1); assert(x0 <= x);
	assert(x1 >= 0); assert(x1 <= 1); assert(x1 >= x);
	float ix = ((x - x0)*c + (x1 - x)*b)/(x1 - x0);
	assert(ix >= a); assert(ix <= c);
	return 1*ix;
}

static
float march_cyclic(float a, float b, float c, float d, float x, float y)
{
	assert(a <= b); assert(b <= d); assert(d <= c);

	// first triangle
	float alpha = c == a ? 0.5 : (b - a)/(c - a);
	if (x <= alpha*(1 - y)) {
		float ix = a;
		if (alpha > 0) {
			float y0=(x+alpha*y)/alpha;
			assert(y0 >= 0); assert(y0 <= 1);
			ix = b*y0 + a*(1 - y0);
		}
		return ix;
	}

	// second triangle
	float beta = c == a ? 0.5 : (d - a)/(c - a);
	if (x >= beta + (1 - beta)*y) {
		float ix = c;
		if (beta < 1) {
			float y0=(1-x+(1-beta)*y)/(1-beta);
			assert(y0 >= 0); assert(y0 <= 1);
			ix = d*y0 + c*(1 - y0);
		}
		return ix;
	}

	// quadrilater
	float x0 = alpha*(1 - y);
	float x1 = beta + (1-beta)*y;
	assert(x0 >= 0); assert(x0 <= 1); assert(x0 <= x);
	assert(x1 >= 0); assert(x1 <= 1); assert(x1 >= x);
	float ix = ((x - x0)*d + (x1 - x)*b)/(x1 - x0);

	return ix;
}


static
float march_singular_raw(float a, float b, float c, float d, float x, float y)
{
	assert(a <= d); assert(d <= b); assert(b <= c);

	// first triangle
	float alpha = c == a ? 0.5 : (b - a)/(c - a);
	if (x <= alpha*(1 - y)) {
		float ix = a;
		if (alpha > 0) {
			float y0=(x+alpha*y)/alpha;
			assert(y0 >= 0); assert(y0 <= 1);
			ix = b*y0 + a*(1 - y0);
		}
		return ix;
	}

	// second triangle
	float gamma = d == c ? 0.5 : (b - c)/(d - c);
	if (y >= 1 - (1 - gamma)*x) {
		float y0 = y+(1-gamma)*(x-1);
		assert(y0 >= 0); assert(y0 <= 1);
		float ix = d*y0 + c*(1 - y0);
		return ix;
	}

	y = 1-y; gamma = 1-gamma;
	if ((y-1)*(1-alpha)<(x-alpha)*(gamma-1)) // constant part
		return b;
	else { // remaining triangle
		float y0 = y+(1-x)*((gamma-1)/(1-alpha));
		y0 = (y0 - gamma)/(1-gamma);
		if (y0 < 0 || y0 >1) return 0;
		assert(y0 >= 0); assert(y0 <= 1);
		float ix = c*y0 + b*(1 - y0);
		return ix;
	}

	return 0;
}

//SMART_PARAMETER(MARCH_SADDLES,2)

static
float march_singular(float a, float b, float c, float d, float x, float y)
{
	int option = 3;//MARCH_SADDLES();
	float mv = (a+b+c+d)/4;
	float d0 = fabs(b - mv);
	float d1 = fabs(d - mv);
	if (0 == option) option = d1<d0 ? -1 : 1;
	if (2 == option) option = d1>d0 ? -1 : 1;
	if (3 == option) option = d1>=d0 ? -1 : 1;
	switch(option) {
	case 1: return march_singular_raw(a,b,c,d,x,y);
	case -1:return -march_singular_raw(-c,-d,-a,-b,1-x,y);
	case -2:return -1;
	default: assert(0);//error("bad MARCH_SADDLES %d", option);
	}
}

static float single_marcho(float a, float b, float c, float d, float x, float y)
{
	assert(x >= 0); assert(y >= 0); assert(x <= 1); assert(y <= 1);
	assert(a <= b); assert(a <= c); assert(a <= d);
	assert(b <= c);

	if (d < b) {
		assert(d <= c);
		return march_singular(a, b, c, d, x, y);
	}
	if (d > c) {
		assert(d >= b);
		return march_regular(a, b, c, d, x, y);
	}
	assert(b <= d);
	assert(d <= c);
	return march_cyclic(a, b, c, d, x, y);
}

float marchi(float a, float b, float c, float d, float x, float y)
{
	if (a <= b && a <= c && a <= d) {
	       if (b <= c)
		       return single_marcho(a, b, c, d, x, y);
	       else // swap diagonal
		       return marchi(a, c, b, d, y, x);
	}
	else if (c <= a && c <= b && c <= d) // swap x axis
		return marchi(c, d, a, b, 1-x, y);
	else if (b <= a && b <= c && b <= d) // swap y axis
		return marchi(b, a, d, c, x, 1-y);
	else if (d <= a && d <= b && d <= c) // swap antidiagonal
		return marchi(d, b, c, a, 1-y, 1-x);
	else
		assert(0);
}
