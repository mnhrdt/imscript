#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

#ifdef __linux
#  include <sys/types.h>
#  include <unistd.h>
static const char *emptystring = "";
static const char *myname(void)
{
#  define n 0x29a
	//const int n = 0x29a;
	static char buf[n];
	pid_t p = getpid();
	snprintf(buf, n, "/proc/%d/cmdline", p);
	FILE *f = fopen(buf, "r");
	if (!f) return emptystring;
	int c, i = 0;
	while ((c = fgetc(f)) != EOF && i < n) {
#  undef n
		buf[i] = c ? c : ' ';
		i += 1;
	}
	if (i) buf[i-1] = '\0';
	fclose(f);
	return buf;
}
#else
static const char *myname(void) { return ""; }
#endif//__linux


static void error(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nERROR(\"%s\"): ", myname());
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#ifdef NDEBUG
	exit(-1);
#else//NDEBUG
	//print_trace(stderr);
	exit(*(int *)0x43);
#endif//NDEBUG
}


typedef float (*extension_operator_float)(float*,int,int,int,int,int,int);

static float extend_float_image_by_zero(float *xx, int w, int h, int pd,
		int i, int j, int l)
{
	float (*x)[w][pd] = (void*)xx;
	if (i < 0 || j < 0 || i > w-1 || j > h-1 || l < 0 || l > pd-1)
		return 0;
	else
		return x[j][i][l];
}

int main(int c, char *v[])
{
	if (c != 4 && c != 5 && c != 6) {
		fprintf(stderr, "usage:\n\t%s w h pd [in [out]]\n", *v);
		return 1;
	}
	int w, ow = atoi(v[1]);
	int h, oh = atoi(v[2]);
	int pd, opd = atoi(v[3]);
	char *fname_in = c > 4 ? v[4] : "-";
	char *fname_out = c > 5 ? v[5] : "-";

	//int w, h, pd;
	//float *x = iio_read_image_float_vec(fname, &w, &h, &pd);
	//printf("%d %d %d\n", w, h, pd);


	float *x = iio_read_image_float_vec(fname_in, &w, &h, &pd);
	float (*y)[ow][opd] = malloc(ow * oh * opd * sizeof(float));
	if (!y) error("out of mem");

	extension_operator_float p = extend_float_image_by_zero;
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < opd; l++)
		y[j][i][l] = p(x, w, h, pd, i, j, l);

	iio_save_image_float_vec(fname_out, y[0][0], ow, oh, opd);

	free(x);
	free(y);

	return 0;
}
