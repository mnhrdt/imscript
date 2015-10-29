#include <assert.h>
#include <stdio.h>
#include "xfopen.c"
#include "parsenumbers.c"

static int getlinen(char *l, int n, FILE *f)
{
	int c, i = 0;
	while (i < n-1 && (c = fgetc(f)) != EOF && c != '\n')
		l[i++] = c;
	l[i] = '\0';
	return i;
}

static void apply_transform(double y[3], double A[12], double x[3])
{
	y[0] = (A[0]*x[0] + A[1]*x[1] + A[2]*x[2] )+ A[3];
	y[1] = (A[4]*x[0] + A[5]*x[1] + A[6]*x[2] )+ A[7];
	y[2] = (A[8]*x[0] + A[9]*x[1] + A[10]*x[2])+ A[11];
}

int main(int c, char *v[])
{
	if (c != 2 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s \"a1 .. a12\" [in [out]]\n", *v);
		//                         0   1            2   3
		return 1;
	}
	char *affstring    = v[1];
	char *filename_in  = c > 2 ? v[2] : "-";
	char *filename_out = c > 3 ? v[3] : "-";

	double A[12];
	read_n_doubles_from_string(A, affstring, 12);
	fprintf(stderr, "3d motion:\n\t%g %g %g %g\n"
			"\t%g %g %g %g\n\t%g %g %g %g\n",
			A[0], A[1], A[2], A[3], A[4], A[5],
			A[6], A[7], A[8], A[9], A[10], A[11]);

	FILE *fi = xfopen(filename_in, "r");
	FILE *fo = xfopen(filename_out, "w");

	int n, nmax = 10000;
	char buf[nmax];
	int cx = 0, maxconvert = 0;
	while (n = getlinen(buf, nmax, fi))
	{
		if (!maxconvert && buf==strstr(buf, "element vertex "))
		{
			int tt, rr = sscanf(buf, "element vertex %d", &tt);
			if (rr == 1)
			{
				maxconvert = tt;
				fprintf(stderr, "treating at most %d points\n",
						maxconvert);
			}
		}
		double x[3], y[3];
		int d, r = sscanf(buf, "%lf %lf %lf %n", x, x+1, x+2, &d);
		if (r == 3 && (!maxconvert || cx < maxconvert)) {
			apply_transform(y, A, x);
			fprintf(fo, "%lf %lf %lf ", y[0], y[1], y[2]);
			fputs(buf + d, fo);
			cx += 1;
		} else
			fputs(buf, fo);
		fputc('\n', fo);
	}
	fprintf(stderr, "converted %d points\n", cx);
	if (maxconvert)
		assert(maxconvert == cx);

	xfclose(fo);
	xfclose(fi);
	return 0;
}
