#define _POSIX_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/mman.h>

#include <X11/Xlib.h>

#include <signal.h>


#include "ftr_x11.c"

static void *alloc_shareable_memory(size_t n)
{
	// note MAP_ANONYMOUS is nicer, but less portable than /dev/zero
#ifdef I_CAN_HAS_MMAP_ANONYMOUS
	int prot = PROT_READ | PROT_WRITE;
	int flags = MAP_SHARED | MAP_ANONYMOUS;
	void *r = mmap(NULL, n, prot, flags, -1, 0);
	if (r == MAP_FAILED) { perror("mmap"); exit(43); }
	return r;
#else
	int fd = open("/dev/zero", O_RDWR);
	if (fd == -1) { perror("open(\"/dev/zero\")"); exit(42); }
	void *r = mmap(NULL, n, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (r == MAP_FAILED) { perror("mmap"); exit(43); }
	if (close(fd)) { perror("close(\"/dev/zero\")"); exit(44); }
	return r;
#endif
}

static struct _FTR *global_variable_for_the_child;

static void my_wait(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	waitpid(f->child_pid, 0, 0);
}

void my_sigusr1_handler(int s)
{
	struct _FTR *f = global_variable_for_the_child;

	fprintf(stderr, "I  got signal %d (%dx%d)\n", s, f->w, f->h);
	fprintf(stderr, "first pixel is %d %d %d\n",
			f->rgb[0], f->rgb[1], f->rgb[2]);

	XEvent ev;
	ev.type = Expose;
	XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
	XFlush(f->display);
}

struct FTR *my_fork_window(int w, int h)
{
	struct FTR *f = alloc_shareable_memory(sizeof*f);

	// place to store the framebuffer for both processes
	void *data = alloc_shareable_memory(3*w*h);

	//*f = ftr_new_window_with_image_uint8_rgb(NULL, w, h);
	//free(f->rgb);
	//f->rgb = alloc_shareable_memory(3*w*h);

	//struct FTR f = ftr_new_window_with_image_uint8_rgb(NULL, w, h);
	//free(f.rgb);
	//f.rgb = alloc_shareable_memory(3*w*h);

	pid_t p = fork();
	if (p) // I'm the parent
	{
		fprintf(stderr,"P: I'm a process with pid %d\n",(int)getpid());
		fprintf(stderr, "P: I've just forked %d\n", (unsigned int)p);
		// todo: stop parent until child starts loop
		//sleep(1);
		int pstatus;
		fprintf(stderr, "P: going to wait for child to send signal\n");
		waitpid(p, &pstatus, WUNTRACED);
		fprintf(stderr, "P: child has stopped (status=%d)\n", pstatus);
		fprintf(stderr, "P: telling the child he can go along\n");
		kill(p, SIGCONT);

		// save the pid of the child
		((struct _FTR*)f)->child_pid = p;
		//f->w = w;
		//f->h = h;
		//f->rgb = data;
		// observation: the parent does NOT have an X connection
		return f;
	}
	else // I'm the child
	{
		fprintf(stderr,"C: I'm a process with pid %d\n",(int)getpid());
		fprintf(stderr, "C: My parent has pid %d\n", (int)getppid());
		*f = ftr_new_window_with_image_uint8_rgb(NULL, w, h);
		free(f->rgb);
		f->rgb = data;

		// set signal handler
		struct sigaction a;
		sigemptyset(&a.sa_mask);
		a.sa_flags = 0;
		a.sa_handler = my_sigusr1_handler;
		sigaction(SIGUSR1, &a, NULL);
		global_variable_for_the_child = (struct _FTR*)f;

		// todo, stop ourselves, so that the parent will
		// know that he can go on
		fprintf(stderr, "C: going to stop myself\n");
		kill(getpid(), SIGSTOP);
		fprintf(stderr, "C: i've somehow unstopped\n");
		exit(ftr_loop_run(f));
	}
}

void my_notify_change(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	fprintf(stderr, "sending signal %d to process %d\n",
			SIGUSR1, (unsigned int)f->child_pid);
	fprintf(stderr, "first pixel is %d %d %d\n",
			f->rgb[0], f->rgb[1], f->rgb[2]);
	kill(f->child_pid, SIGUSR1);

//	ev.type = Expose;
//	XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
//	XFlush(f->display);
}

int main_forks(int c, char *v[])
{
	struct FTR *f = my_fork_window(320,200);
	struct FTR *g = my_fork_window(320,200);
	for (int i = 0; i < f->w * f->h; i++) {
		f->rgb[3*i+0] = 0; f->rgb[3*i+1] = 255; f->rgb[3*i+2] = 0;
		g->rgb[3*i+0] = 0; g->rgb[3*i+1] = 255; g->rgb[3*i+2] = 127;
	}
	my_notify_change(f);
	my_notify_change(g);
	sleep(1);
	for (int i = 0; i < f->w * f->h; i++) {
		f->rgb[3*i+0] = 255; f->rgb[3*i+1] = 0; f->rgb[3*i+2] = 0;
		g->rgb[3*i+0] = 255; g->rgb[3*i+1] = 0; g->rgb[3*i+2] = 127;
	}
	my_notify_change(f);
	my_notify_change(g);
	sleep(1);
	for (int i = 0; i < f->w * f->h; i++) {
		f->rgb[3*i+0] = 0; f->rgb[3*i+1] = 0; f->rgb[3*i+2] = 255;
		g->rgb[3*i+0] = 0; g->rgb[3*i+1] = 0; g->rgb[3*i+2] = 127;
	}
	my_notify_change(f);
	my_notify_change(g);
	sleep(1);
	for (int i = 0; i < f->w * f->h; i++) {
		f->rgb[3*i+0] = 255; f->rgb[3*i+1] = 255; f->rgb[3*i+2] = 255;
		g->rgb[3*i+0] = 255; g->rgb[3*i+1] = 255; g->rgb[3*i+2] = 127;
	}
	my_notify_change(f);
	my_notify_change(g);
	my_wait(f);
	my_wait(g);
	return 0;
}

int main(int c, char *v[]) { return main_forks(c, v); }
