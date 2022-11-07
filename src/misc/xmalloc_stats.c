#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <execinfo.h>

#define BACKTRACE_SYMBOLS 50

void print_trace(FILE *f)
{
	void *array[BACKTRACE_SYMBOLS];
	size_t size, i;
	char **strings;

	size = backtrace (array, BACKTRACE_SYMBOLS);
	strings = backtrace_symbols (array, size);

	fprintf (f, "Obtained %zu stack frames.\n", size);

	for (i = 0; i < size; i++)
		fprintf (f, "%s\n", strings[i]);

	free (strings);
}

void bad_fail(const char *fmt, ...)
{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
	print_trace(stderr);
	exit(*(volatile int *)0);
}


#define SIZE_OF_MEGABYTE 1048576.0
static size_t xmalloc_global_size = 0;
static size_t xmalloc_global_id = 0;
struct xmalloc_secret_stuff { size_t size; int id; };
static struct xmalloc_secret_stuff *xmalloc_get_secret_stuff(void *x)
{
	return (void *)(-sizeof(struct xmalloc_secret_stuff) + (char *)x);
}
static void *xmalloc_get_stuff(struct xmalloc_secret_stuff *p)
{
	return sizeof*p + (char *)p;
}
void xmalloc_debug(const char *fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
}


void print_malloc_stats(void)
{
	xmalloc_debug(" xgs = %zu (%gMB)\n", xmalloc_global_size,
			xmalloc_global_size/SIZE_OF_MEGABYTE);
}



static void xmalloc_error(size_t s)
{
	double sm = s / (0x100000 * 1.0);
	bad_fail("xmalloc: out of memory while requesting "
			"%zu bytes (%gMB):\"%s\"",
			s, sm, strerror(errno));
}

void *xmalloc(size_t size)
{
	void *new;

	size_t id = xmalloc_global_id;
	xmalloc_debug("XMALLOC:\tid = %.3zu, size = %8zu (%10.5gMb) ", id, size,
			size/SIZE_OF_MEGABYTE);

	if (size == 0)
		bad_fail("xmalloc: zero size");
	//new = malloc(size + 2 * sizeof(size_t));
	new = malloc(size + sizeof(struct xmalloc_secret_stuff));
	if (new == NULL)
		xmalloc_error(size);
		//error("xmalloc: out of memory while requesting "
		//		"%d bytes (%s)",
		//		size, strerror(errno));

	xmalloc_debug("{%p}\t\t\t\t\t\t", new);

	struct xmalloc_secret_stuff *secret = new;
	void *stuff = xmalloc_get_stuff(secret);
	assert(secret == xmalloc_get_secret_stuff(stuff));

	xmalloc_global_size += size;
	xmalloc_global_id += 1;
	secret->id = id;
	secret->size = size;
	//*(size_t *)new = id;
	//*(1+(size_t *)new) = size;
	print_malloc_stats();
	//xmalloc_debug("XMALLOC2:\t\tid= %.3zu\n", id);

	//return (char *)new + 2 * sizeof(size_t);
	return stuff;
}

void xfree(void *p)
{
	if (!p)
		bad_fail("thou shalt not free a null pointer!");
	struct xmalloc_secret_stuff *secret = xmalloc_get_secret_stuff(p);
	assert(p == xmalloc_get_stuff(secret));
	xmalloc_global_size -= secret->size;
	xmalloc_debug("XFREE:\t\tid = %.3zu, size = %8zu (%10.5gMb)\t\t\t\t\t\t\t", secret->id, secret->size, secret->size/SIZE_OF_MEGABYTE);

	print_malloc_stats();

	free(secret);
}

