#include <stdbool.h>

//
// insert figure "fours" here
//


static struct {
	int id;
	int v[2];
	int d[2];
	int e[2];
	//int neighbor[2];
} table_of_fours[] = {
	{ 0, {0,0},  {1,0}, {1,3} },
	{ 1, {1,0},  {0,1}, {3,2} },
	{ 2, {0,1}, {-1,0}, {2,0} },
	{ 3, {1,1}, {0,-1}, {0,1} },
	{ -1, {0,0}, {0,0}, {0,0} }
};

static bool inbetweenP(float x, float y, float m)
{
	if (x <= m && m <= y) return true;
	if (y <= m && m <= x) return true;
	return false;
}

static int march_one_step_from_below(float v[4], float t)
{
	assert(inbetweenP(v[0], v[1], t));
	int k = marching_squares_case(v[0], v[1], v[2], v[3]);
	assert(marching_squares_table[k].frombelow);
	return marching_squares_table[k].belowto;
}

// when entering a square from edge indexed "from", where do we exit?
int march_one_step(float values[4], float threshold, int from)
{
	assert(from >= 0 && from < 4);
	int k = marching_squares_case(a, b, c, d, t);
	if (k == 0 || k == 15) error("empty case!");
	float w[4];
	FORI(4)
		w[i] = v[(i-from)%4];
	int r = march_one_step_from_below(w, threshold);
	return (r+from)%4;
}
