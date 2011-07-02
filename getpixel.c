
typedef float (*getsample_operator)(float*,int,int,int,int,int,int);
//typedef void (*setsample_operator)(float*,int,int,int,int,int,int,float);

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

static float getsample_error(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return *(int*)0;
	return x[(i+j*w)*pd + l];
}

static void setsample_0(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
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
