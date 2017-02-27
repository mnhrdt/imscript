#ifdef EXTRAPOLATORS_C
#define EXTRAPOLATORS_C
// A "extrapolator" evaluates an image at an arbitrary integral position.
// When the position is outside the image domain, the value is extrapolated
// (by periodization, or by a constant value).

// type of the "extrapolator" functions
typedef float (*extrapolator_t)(float*,int,int,int,int,int,int);

// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int n, int p)
{
	int r = n % p;
	r = r < 0 ? r + p : r;
	assert(r >= 0);
	assert(r < p);
	return r;
}

// instance of "extrapolator_t", extrapolate by periodicity
static float getpixel_per(float *x, int w, int h, int i, int j)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	return x[i+j*w];
}

// instance of "extrapolator_t", extrapolate by a constant value
static float getpixel_cons(float *x, int w, int h, int i, int j)
{
	static float value = 0;
	if (w == 0 && h == 0) // if the image has zero size, set the constant
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return value;
	return x[i+j*w];
}

// instance of "extrapolator_t", extrapolate by a nearest neighbor
static float getpixel_nn(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[i+j*w];
}

// auxiliary function: 2p-periodic replication
static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p) r -= 1;
	if (r > p)
		r = 2*p - r;
	if (n < 0 && p > 1) r += 1;
	assert(r >= 0);
	assert(r < p);
	return r;
}

// instance of "extrapolator_t", extrapolate by a nearest neighbor
static float getpixel_sym(float *x, int w, int h, int i, int j)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	return x[i+j*w];
}

#endif//EXTRAPOLATORS_C
