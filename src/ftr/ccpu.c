// CCPU : a controller for cpu


#include <assert.h>
#include <signal.h>  // pid_t, kill
#include <stdlib.h>  // system
#include <stdio.h>   // snprintf, fprintf
#include <spawn.h>   // posix_spawn
#include <unistd.h>  // getpid
#include "iio.h"     // iio_image_write

#define MAX_CPUS 10

// internal struct
struct cpu_view {
	pid_t p;  // process id of this cpu viewer
	char f[FILENAME_MAX];   // interchange filename
};

// internal data
struct cpu_view global_table_of_cpu_views[MAX_CPUS];


#define CCPU_SHOW_DEBUG_MESSAGES

#ifdef CCPU_SHOW_DEBUG_MESSAGES
#  define CCPU_DEBUG(...) do {\
	if (getenv("CCPU_DEBUG")) fprintf(stderr, __VA_ARGS__); } while(0)
#else//CCPU_SHOW_DEBUG_MESSAGES
#  define CCPU_DEBUG(...) do;while(0) /* nothing */
#endif// CCPU_SHOW_DEBUG_MESSAGES


// API: create a new view with the given image buffer
// returns a handle to the view
int cpu_new(float *x, int w, int h, int d)
{
	struct cpu_view *v = global_table_of_cpu_views + 0; // TODO: do increase this

	snprintf(v->f, FILENAME_MAX,
			"/tmp/cpu_view_%d_%d.npy", getuid(), getpid());
	CCPU_DEBUG("cpu_new f = %s\n", v->f);
	iio_write_image_float_vec(v->f, x, w, h, d);
	char *a[] = {"cpu", v->f, NULL};
	extern char **environ;
	int r = posix_spawnp(&v->p, a[0], NULL, NULL, a, environ);
	CCPU_DEBUG("cpu_new r = %d\n", r);
	CCPU_DEBUG("cpu_new p = %d\n", v->p);

	return 0;
}

// API: send a key to a cpu
void cpu_send_key(int n, int k)
{
	CCPU_DEBUG("cpu_send_key n=%d k=%d\n", n, k);
	assert(n == 0);
	struct cpu_view *v = global_table_of_cpu_views + n;

	char c[FILENAME_MAX]; // command line to run
	snprintf(c, FILENAME_MAX,
		"xdotool search --any --pid %d --name ftr_win_pid_%d key %c",
		v->p, v->p, k);
	CCPU_DEBUG("cpu_send_key c=\"%s\"\n", c);
	// NOTE: this xdotool call is a subtle hack to support transparently
	// plain x11 and glut FTR windows.  The "pid" field is only used for
	// plain x11, the "name" field is only used for ftr_freeglut.
	int r = system(c);
	(void)r;
}

// API: update the view determined by this handle
void cpu_update(int n, float *x, int w, int h, int d)
{
	CCPU_DEBUG("cpu_update n=%d x=%p %dx%d,%d\n", n, (void*)x, w, h, d);
	assert(n == 0);
	struct cpu_view *v = global_table_of_cpu_views + n;

	if (x)
		iio_write_image_float_vec(v->f, x, w, h, d);
	cpu_send_key(n, '2');
	// TODO: maybe send an actual expose event to the window?
	// caveat: i don't know how to capture it from glut, so it would need
	// some more hacking on the client side
}

// API: close the cpu window
void cpu_close(int n)
{
	CCPU_DEBUG("cpu_close n=%d\n", n);
	assert(n == 0);
	struct cpu_view *v = global_table_of_cpu_views + n;

	cpu_send_key(n, 'q');
	unlink(v->f);
}

int main_try_ccpu(void)
{
	int w, h, d;
	float *x = iio_read_image_float_vec("/tmp/barbara.png", &w, &h, &d);
	int n = cpu_new(x, w, h, d);
	free(x);
	sleep(5);
	x = iio_read_image_float_vec("/tmp/lenak.png", &w, &h, &d);
	cpu_update(n, x, w, h, d);
	sleep(3);
	cpu_close(n);
	return 0;
}

#ifdef MAIN_CCPU
int main(void)
{
	return main_try_ccpu();
}
#endif
