#include <assert.h>

#ifndef xmalloc
#include <stdlib.h>
#define xmalloc malloc
#endif//xmalloc


#define V00 0 //{0,0}
#define V10 1 //{1,0}
#define V01 2 //{0,1}
#define V11 3 //{1,1}
#define VNO -1 //{-1,-1}
#define E0 {V00,V10}
#define E1 {V10,V11}
#define E2 {V11,V01}
#define E3 {V01,V00}
#define E_ {VNO,VNO}
static const struct { int c,ns; int s[2][2][2]; } marching_squares_table[] = {
	{ 0, 0, {{E_,E_},{E_,E_}} },
	{ 1, 1, {{E0,E3},{E_,E_}} },
	{ 2, 1, {{E0,E1},{E_,E_}} },
	{ 3, 1, {{E1,E3},{E_,E_}} },
	{ 4, 1, {{E3,E2},{E_,E_}} },
	{ 5, 1, {{E0,E2},{E_,E_}} },
	{ 6, 2, {{E1,E2},{E3,E0}} },
	{ 7, 1, {{E1,E2},{E_,E_}} },
	{ 8, 1, {{E2,E1},{E_,E_}} },
	{ 9, 2, {{E0,E1},{E2,E3}} },
	{10, 1, {{E2,E0},{E_,E_}} },
	{11, 1, {{E2,E3},{E_,E_}} },
	{12, 1, {{E3,E1},{E_,E_}} },
	{13, 1, {{E0,E1},{E_,E_}} },
	{14, 1, {{E3,E0},{E_,E_}} },
	{15, 0, {{E_,E_},{E_,E_}} }
};

static float vertex_table[5][2] = { {0,0}, {1,0}, {0,1}, {1,1}, {-1,-1} };

static void interpolate_two_points(float out[2], float a[2], float b[2],
		float va, float vb, float t)
{
	assert(va != vb);
	for (int i = 0; i < 2; i++)
		out[i] = ((vb - t)*a[i] + (t - va)*b[i])/(vb - va);
}

/*
 *      v[2] ---- v[3]
 *       |         |
 *       |         |
 *       |         |
 *      v[0] ---- v[1]
 *
 */
int marching_squares_cell(float (*segments)[2][2], float v[4], float t)
{
	int m = 0;
	for (int i = 0; i < 4; i++)
		if (v[i] >= t)
			m += (1 << i);
	int r = marching_squares_table[m].ns;
	for (int i = 0; i < r; i++) {
		for (int k = 0; k < 2; k++) {
			int p[2];
			float *f[2], g[2];
			for (int j = 0; j < 2; j++) {
				p[j] = marching_squares_table[m].s[i][k][j];
				f[j] = vertex_table[p[j]];
			}
			interpolate_two_points(g,f[0],f[1],v[p[0]],v[p[1]],t);
			for (int j = 0; j < 2; j++)
				segments[i][k][j] = g[j];
		}
	}
	return r;
}


float (*marching_squares_whole_image_float(int *n,
			float *image, int w, int h, float t))[2][2]
{
#define getpixel(x,w,i,j) (x)[(i)+(w)*(j)]
	static const int V[4][2] = { {0,0}, {1,0}, {0,1}, {1,1} };
	float (*segments)[2][2] = malloc(w*h*2*sizeof*segments);
	*n = 0;
	for (int j = 0; j < h-1; j++)
		for (int i = 0; i < w-1; i++) {
			float ts[2][2][2], v[4];
			for (int l = 0; l < 4; l++)
				v[l] = getpixel(image,w, i+V[l][0], j+V[l][1]);
			int ns = marching_squares_cell(ts, v, t);
			for (int k = 0; k < ns; k++) {
				segments[*n][0][0] = ts[k][0][0] + i;
				segments[*n][0][1] = ts[k][0][1] + j;
				segments[*n][1][0] = ts[k][1][0] + i;
				segments[*n][1][1] = ts[k][1][1] + j;
				*n += 1;
			}
		}
	return segments;
#undef getpixel
}



#ifdef USE_MAIN

#include <stdio.h>

int main(int c, char *v[])
{
	int n;
	float x[9] = {0,0,0, 0,255,0, 0,0,0};
	float (*s)[2][2] = marching_squares_whole_image_float(&n, x, 3, 3,
			127.5);
	for (int i = 0; i < n; i++)
		printf("segment %d of %d: (%g %g) - (%g %g)\n", i, n,
				s[i][0][0], s[i][0][1],
				s[i][1][0], s[i][1][1]);
	return 0;
}
#endif//USE_MAIN
