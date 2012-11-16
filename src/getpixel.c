#ifndef _GETPIXEL_C
#define _GETPIXEL_C

typedef float (*getsample_operator)(float*,int,int,int,int,int,int);
//typedef void (*setsample_operator)(float*,int,int,int,int,int,int,float);

// extrapolate by 0
inline
static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

// extrapolate by nearest value
inline
static float getsample_1(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	return x[(i+j*w)*pd + l];
}

inline
static float getsample_error(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return *(volatile int*)0;
	return x[(i+j*w)*pd + l];
}

inline
static void setsample_0(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

inline
static void setsample_segf(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		*(volatile int*)0 = 0;
	x[(i+j*w)*pd + l] = v;
}

//static float getpixel(float *x, int w, int h, int i, int j)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h)
//		return 0;
//	return x[i + j*w];
//}
//
//static void setpixel(float *x, int w, int h, int i, int j, float v)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h)
//		return;
//	x[i + j*w] = v;
//}

#endif//_GETPIXEL_C
