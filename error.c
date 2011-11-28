#ifndef _ERROR_C
#define _ERROR_C

#include <stdio.h>

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


#include <stdarg.h>

static void error(const char *fmt, ...) __attribute__((noreturn,format(printf,1,2)));
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

#endif//_ERROR_C
