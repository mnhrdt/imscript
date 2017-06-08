#include <math.h>
#include <stdio.h>

float holdwater(float *H0, int n)
{
	// preserve original data
	float H[n];
	for (int i = 0; i < n; i++)
		H[i] = H0[i];

	// find summit
	float max = -INFINITY;
	int imax = 0;
	for (int i = 0; i < n; i++)
		if (H[i] > max) {
			max = H[i];
			imax = i;
		}

	// fill before summit
	int lasttop = 0;
	for (int i = 0; i < imax; i++)
		if (H[lasttop] > H[i])
			H[i] = H[lasttop];
		else
			lasttop = i;

	// fill after summit
	lasttop = n - 1;
	for (int i = n - 1; i > imax; i--)
		if (H[lasttop] > H[i])
			H[i] = H[lasttop];
		else
			lasttop = i;

	// compute difference between the filled and unfilled versions
	float total = 0;
	for (int i = 0; i < n; i++)
		total += H[i] - H0[i];

	return total;
}

float hw(float *H, int n)
{
	float r = -INFINITY, t = *H;
	int i, j = 0, k = 0;

	for (i = 0; i < n; i++)
		if (H[i] > r)
			r = H[k=i];

	r = j = 0;
	for (i = 0; i < k; i++)
		if (t > H[i])
			r += t - H[i];
		else
			t = H[j=i];

	t = H[j=n-1];
	for (i = n - 1; i > k; i--)
		if (t > H[i])
			r += t - H[i];
		else
			t = H[j=i];

	return r;
}

int main(void)
{
	float H[] = {1,3,1,3,1,2,0,1,3,1,2,1,1};
	float w = holdwater(H, sizeof H / sizeof *H);
	float w2 = hw(H, sizeof H / sizeof *H);
	printf("%g\n", w);
	printf("%g\n", w2);
	return 0;
}
