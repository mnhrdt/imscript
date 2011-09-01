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


static char *typstring(int w, bool ieeefp, bool sign)
{
	if (ieeefp) {
		if (w == 4) return "float";
		if (w == 8) return "double";
		error("uncecognized ieee floating point width %d", w);
	} else {
		if (w == 1) return sign ? "int8" : "uint8";
		if (w == 2) return sign ? "int16" : "uint16";
		if (w == 4) return sign ? "int32" : "uint32";
		if (w == 8) return sign ? "int64" : "uint64";
		error("uncecognized integer width %d", w);
	}
	return "x";
}

int main(int c, char *v[])
{
	if (c != 1 && c != 2) {
		fprintf(stderr, "usage:\n\t%s [in]\n", *v);
		return 1;
	}
	char *fname = c > 1 ? v[1] : "-";

	//int w, h, pd;
	//float *x = iio_read_image_float_vec(fname, &w, &h, &pd);
	//printf("%d %d %d\n", w, h, pd);

	int dim, siz[IIO_MAX_DIMENSION], pd, ss;
	bool ieeefp, sgn;
	void *x = iio_read_nd_image_as_stored(fname, &dim, siz, &pd, &ss, &ieeefp, &sgn);

	if (!x) { printf("unreadable\n"); return 2; }
	char *s = typstring(ss, ieeefp, sgn);
	if (dim == 2) {
		printf("%d %d %d %s\n", siz[0], siz[1], pd, s);
	} else {
		printf("%d", dim);
		for (int i = 0; i < dim; i++)
			printf(" %d", siz[i]);
		printf(" %d %s\n", pd, s);
	}

	return 0;
}
