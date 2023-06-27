// CCPU : a controller for cpu


#include <assert.h>
#include <signal.h>  // pid_t, kill
#include <stdlib.h>  // system
#include <stdio.h>   // snprintf, fprintf
#include <spawn.h>   // posix_spawn
#include <unistd.h>  // getpid
#include "iio.h"     // iio_image_write

#define MAX_CPUS 10

struct cpu_view {
	pid_t p;  // process id of this cpu viewer
	char f[FILENAME_MAX];   // interchange filename
};

struct cpu_view global_table_of_cpu_views[MAX_CPUS];


// API: create a new view with the given image buffer
// returns a handle to the view
int cpu_new(float *x, int w, int h, int d)
{
	struct cpu_view *v = global_table_of_cpu_views + 0;

	snprintf(v->f, FILENAME_MAX,
			"/tmp/cpu_view_%d_%d.npy", getuid(), getpid());
	fprintf(stderr, "cpu_new f = %s\n", v->f);
	iio_write_image_float_vec(v->f, x, w, h, d);
	char *a[] = {"cpu", v->f, NULL};
	extern char **environ;
	int r = posix_spawnp(&v->p, a[0], NULL, NULL, a, environ);
	fprintf(stderr, "cpu_new r = %d\n", r);
	fprintf(stderr, "cpu_new p = %d\n", v->p);

	return 0;
}

void cpu_send_key(int n, int k)
{
	assert(n == 0);
	struct cpu_view *v = global_table_of_cpu_views + n;

	char c[FILENAME_MAX]; // command line to run
	snprintf(c, FILENAME_MAX,
		"xdotool search --any --pid %d --name ftr_win_pid_%d key %c",
		v->p, v->p, k);
	system(c);
}

// API: update the view determined by this handle
void cpu_update(int n, float *x, int w, int h, int d)
{
	assert(n == 0);
	struct cpu_view *v = global_table_of_cpu_views + n;

	iio_write_image_float_vec(v->f, x, w, h, d);
	cpu_send_key(n, '2');
}

// API: close the cpu window
void cpu_close(int n)
{
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
