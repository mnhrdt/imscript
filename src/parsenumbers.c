#ifndef _PARSENUMBERS_C
#define _PARSENUMBERS_C

// utility function: parse floats from a text file
// returns a pointer to a malloc'ed array of the parsed floats
// fills *no with the number of floats
static float *read_ascii_floats(FILE *f, int *no)
{
	int r, n = 0, nt = 0;
	float *t = NULL;
	while(1) {
		if (n >= nt)
		{
			nt = 2 * (nt + 1);
			t = xrealloc(t, nt * sizeof * t);
		}
		r = fscanf(f, "%f", t+n);
		if (r != 1)
			break;
		n += 1;
	}
	*no = n;
	return t;
}

#endif//_PARSENUMBERS_C
