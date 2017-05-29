#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif
#include <unistd.h>
#include <time.h>
double seconds(void)
{
	static int initialized = 0;
	static time_t first_seconds;
	struct timespec t[1];
	if (!initialized) {
		clock_gettime(CLOCK_REALTIME, t);
		first_seconds = t->tv_sec;
		initialized = 1;
	}
	clock_gettime(CLOCK_REALTIME, t);
	//assert(t->tv_sec >= first_seconds);
	double r = (t->tv_sec - first_seconds) + 1e-9*t->tv_nsec;
	return r;
}
