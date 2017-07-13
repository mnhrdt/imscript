#ifndef DRAWTRIANGLE_C
#define DRAWTRIANGLE_C

#include <assert.h>
#include <stdlib.h> // qsort
#include <math.h>   // fmax, fmin, floor, ceil

// cut the line between points a and b at height h
// note: the line must not be horizontal!
static float cut_line_at_height(float a[2], float b[2], float h)
{
	assert(a[1] != b[1]);
	assert((a[1] - h) * (b[1] - h) < 0);
	float t = (h - a[1]) / (b[1] - a[1]);
	float r = (1 - t) * a[0] + t * b[0];
	return r;
}

// used for sorting points according to height
static int compare_points_by_height(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (a[1] < b[1]) - (a[1] > b[1]);
}

// split a triangle into two triangles with and horizontal side
static void split_triangle_horizontally(
		float pqr[3][2],   // flat bottom half
		float xyz[3][2],   // flat upper half
		float abc[3][2]    // input triangle in arbitrary position
		)
{
	qsort(abc, 3, sizeof*abc, compare_points_by_height);
	pqr[0][0] = xyz[0][0] = abc[1][0];
	pqr[0][1] = xyz[0][1] = abc[1][1];
	pqr[1][0] = xyz[1][0] = cut_line_at_height(abc[0], abc[2], abc[1][1]);
	pqr[1][1] = xyz[1][1] = abc[1][1];
	pqr[2][0] = abc[0][0]; pqr[2][1] = abc[0][1];
	xyz[2][0] = abc[2][0]; xyz[2][1] = abc[2][1];
}

// traverse a triangle whose first side is horizontal
static void traverse_flat_triangle(
		float abc[3][2],          // coordinates of the three points
		void (*f)(int,int,void*), // function to apply to each pixel
		void *e                   // passed to f
		)
{
	assert(abc[0][1] == abc[1][1]);
	int j_bot = ceil (fmin(abc[0][1], abc[2][1]));
	int j_top = floor(fmax(abc[0][1], abc[2][1]));
	for (int j = j_bot; j <= j_top; j++)
	{
		float x = cut_line_at_height(abc[0], abc[2], j);
		float y = cut_line_at_height(abc[1], abc[2], j);
		int i_lef = ceil (fmin(x, y));
		int i_rig = floor(fmax(x, y));
		for (int i = i_lef; i <= i_rig; i++)
			f(i, j, e);
	}
}

// fill a triangle defined by three points a,b,c
void traverse_triangle(
		float abc[3][2],          // coordinates of the three points
		void (*f)(int,int,void*), // function to apply to each pixel
		void *e                   // passed to f
		)
{
	float pqr[3][2]; // the flat bottomed half of triangle abc
	float xyz[3][2]; // the upper-flat half of the triangle abc
	split_triangle_horizontally(pqr, xyz, abc);
	traverse_flat_triangle(pqr, f, e);
	traverse_flat_triangle(xyz, f, e);
}

#endif//DRAWTRIANGLE_C
