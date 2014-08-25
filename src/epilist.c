struct epipolar_line {
	float origin[2];
	float direction[2];
	float length;
};

void epiget(struct epipolar_line *e,
		float x[2], long double fm[9], float w, float h)
{
	long double abc[3] = {
		fm[0] * x[0] + fm[1] * x[1] + fm[2],
		fm[3] * x[0] + fm[4] * x[1] + fm[5],
		fm[6] * x[0] + fm[7] * x[1] + fm[8]
	};
	long double nabc = hypot(hypot(abc[0], abc[1]), abc[2]);
	if (abc[1] < 0)
		nabc *= -1;
	abc[0] /= nabc;
	abc[1] /= nabc;
	abc[2] /= nabc;
}

void fill_epilist(struct epipolar_line *e, int w, int h, long double fm[9])
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float x[2] = {i, j};
		struct epipolar_line *ex = e + j*w + i;
		epiget(ex, x, fm, w, h);
	}
}
