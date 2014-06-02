#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h> // only for "fork"
#include <GL/freeglut.h>

#include "ftr.h"

struct _FTR {
	// visible state
	int w, h, max_w, max_h;
	unsigned char *rgb;
	int stop_loop;
	int changed;
	void *userdata;


	// user-supplied handlers (not visible, but common to all backends)
	ftr_event_handler_t handle_key;
	ftr_event_handler_t handle_button;
	ftr_event_handler_t handle_motion;
	ftr_event_handler_t handle_expose;
	ftr_event_handler_t handle_expose2;
	ftr_event_handler_t handle_resize;
	ftr_event_handler_t handle_idle;
	ftr_event_handler_t handle_idle_toggled;

	// glut-only internal data
	// (there's none, glut requires using global variables)
	//int glut_stopped;
	int handle_mute;
	int glut_initted;
	int glut_right_button_pressed;
	
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];


// glut-specific part {{{1

// global variable
struct _FTR *ftr_freeglut_global_state = NULL;

// glut handlers {{{2
static void my_displayfunc(void)
{
	//fprintf(stderr, "GLUT displayfunc\n");

	struct _FTR *f = ftr_freeglut_global_state;

	if (f->handle_expose)
		f->handle_expose((void*)f, 0, 0, 0, 0); 

	glViewport(0, 0, f->w, f->h);

	glRasterPos2i(-1,1);
	glPixelZoom(1,-1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // ugly

	//fprintf(stderr, "glDrawPixels %p\n", f->rgb);
	glDrawPixels(f->w, f->h, GL_RGB, GL_UNSIGNED_BYTE, f->rgb);

	glutSwapBuffers();
}

static void my_reshapefunc(int w, int h)
{
	fprintf(stderr, "GLUT reshapefunc %d %d\n", w, h);
}

static void my_keyboardfunc(unsigned char k, int x, int y)
{
	fprintf(stderr, "GLUT keyboardfunc  %d '%c' %d %d\n",
			k, isalnum(k)?k:' ', x, y);

	struct _FTR *f = ftr_freeglut_global_state;

	if (f->handle_key) {
		f->handle_key((void*)f, k, 0, x, y); 
		if (f->changed)
			glutPostRedisplay();
	}
}

static void my_specialfunc(int k, int x, int y)
{
	fprintf(stderr, "GLUT specialfunc  %d '%c' %d %d\n",
			k, isalnum(k)?k:' ', x, y);

	//if (f->handle_key)
	//	f->handle_key(f, from_glut_key_to_ftr(k), 0, x, y); 
}

static void my_mousefunc(int b, int s, int x, int y)
{
	fprintf(stderr, "GLUT mousefunc %d %d (%d %d)\n", b, s, x, y);

	struct _FTR *f = ftr_freeglut_global_state;

	if (f->handle_button && 0 == s) {
		f->handle_button((void*)f, b+1, 0, x, y); 
		if (f->changed)
			glutPostRedisplay();
	}
}

static void my_motionfunc(int x, int y)
{
	fprintf(stderr, "GLUT motionfunc %d %d\n", x, y);
}

static void my_passivemotionfunc(int x, int y)
{
	struct _FTR *f = ftr_freeglut_global_state;
	if (f->handle_mute) return;

	fprintf(stderr, "GLUT passive motionfunc %d %d\n", x, y);

	f->handle_mute = 1;
	if (f->handle_motion) {
		f->handle_motion((void*)f, 0, 0, x, y); 
		if (f->changed)
			glutPostRedisplay();
	}
	f->handle_mute = 0;
}

static void my_idle(void)
{
	//fprintf(stderr, "GLUT idle handler\n");

	struct _FTR *f = ftr_freeglut_global_state;
	f->handle_idle((void*)f, 0, 0, 0, 0);
	if (f->changed)
		glutPostRedisplay();
}

// setup_glut_environment {{{2
static void setup_glut_environment(struct _FTR *f)
{
	int argc = 1;
	char *argv[]={"dummy"};
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
			GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutCreateWindow("Whatever");
	glutReshapeWindow(f->w, f->h);
	glutDisplayFunc(my_displayfunc);
	glutReshapeFunc(my_reshapefunc);
	glutIdleFunc(NULL);
	glutMouseFunc(my_mousefunc);
	glutMotionFunc(my_motionfunc);
	glutPassiveMotionFunc(my_passivemotionfunc);
	glutKeyboardFunc(my_keyboardfunc);
	glutSpecialFunc(my_specialfunc);
	f->glut_initted = 1;
	fprintf(stderr, "setup glut rgb = %p\n", f->rgb);
}


// ftr_new_window_with_image_uint8_rgb {{{2
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];

	f->w = w;
	f->h = h;
	f->max_w = 2000;
	f->max_h = 2000;
	f->rgb = malloc(f->max_w * f->max_h * 3);
	for (int i = 0; i < 3*w*h; i++)
		f->rgb[i] = x? x[i] : 0;


	f->handle_key = ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->stop_loop = 0;
	f->handle_mute = 0;
	f->glut_initted = 0;

	//f->glut_stopped = 0;
	// freeglut specific part
	setup_glut_environment(f);

	//f->handle_expose2 = ftr_handler_stop_loop;
	//ftr_loop_run((struct FTR *)f);
	//f->handle_expose2 = 0;

	return *(struct FTR *)f;
}

void ftr_close(struct FTR *ff)
{
	// do nothing by now

	//struct _FTR *f = (void*)ff;
	//if (f->ximage) XDestroyImage(f->ximage);
	//if (f->rgb) free(f->rgb);
	//XCloseDisplay(f->display);
}

int ftr_loop_run(struct FTR *ff)
{
	fprintf(stderr, "going to start loop\n");
	struct _FTR *f = (void*)ff;

	// after window creation, we jumped into the air
	// here, we get hold of a new branch
	ftr_freeglut_global_state = f; 

	if (!f->glut_initted)
		setup_glut_environment(f);
	
	if (f->handle_idle)
		glutIdleFunc(my_idle);

	glutMainLoop();
	fprintf(stderr, "returned from glut loop (%d)\n", f->stop_loop);

	int r = f->stop_loop;
	f->stop_loop = 0;
	return r;
}

void ftr_notify_the_desire_to_stop_this_loop(struct FTR *ff, int retval)
{
	fprintf(stderr, "stop notification retval = %d\n", retval);
	struct _FTR *f = (void*)ff;
	f->stop_loop = retval;
	f->glut_initted = 0;
	glutLeaveMainLoop();
}

// common to all implementations {{{1

struct FTR ftr_new_window(int w, int h)
{
	return ftr_new_window_with_image_uint8_rgb(NULL, w, h);
}

void ftr_fork_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	if (!fork())
	{
		struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
		exit(ftr_loop_run(&f));
	}
}

void ftr_loop_fork(struct FTR *f)
{
	if (!fork())
		exit(ftr_loop_run(f));
}

void ftr_handler_exit_on_ESC(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_exit_on_ESC_or_q(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_stop_loop(struct FTR *f, int k, int m, int x, int y)
{
	ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_dummy(struct FTR *f, int k, int m, int x, int y)
{
}

void ftr_handler_toggle_idle(struct FTR *ff, int k, int m, int x, int y)
{
	struct _FTR *f = (void*)ff;
	ftr_event_handler_t tmp = f->handle_idle;
	f->handle_idle = f->handle_idle_toggled;
	f->handle_idle_toggled = tmp;
}

int ftr_set_handler(struct FTR *ff, char *id, ftr_event_handler_t e)
{
	struct _FTR *f = (void*)ff;
	if (0) ;
	else if (0 == strcmp(id, "key"))    { f->handle_key    = e; return 0; }
	else if (0 == strcmp(id, "button")) { f->handle_button = e; return 0; }
	else if (0 == strcmp(id, "motion")) { f->handle_motion = e; return 0; }
	else if (0 == strcmp(id, "expose")) { f->handle_expose = e; return 0; }
	else if (0 == strcmp(id, "resize")) { f->handle_resize = e; return 0; }
	else if (0 == strcmp(id, "idle"))   { f->handle_idle   = e; return 0; }
	return fprintf(stderr, "WARNING: unrecognized event \"%s\"\n", id);
}

ftr_event_handler_t ftr_get_handler(struct FTR *ff, char *id)
{
	struct _FTR *f = (void*)ff;
	if (0) ;
	else if (0 == strcmp(id, "key"))    return f->handle_key   ;
	else if (0 == strcmp(id, "button")) return f->handle_button;
	else if (0 == strcmp(id, "motion")) return f->handle_motion;
	else if (0 == strcmp(id, "expose")) return f->handle_expose;
	else if (0 == strcmp(id, "resize")) return f->handle_resize;
	else if (0 == strcmp(id, "idle"))   return f->handle_idle  ;
	return fprintf(stderr, "WARNING: bad event \"%s\"\n", id),
		(ftr_event_handler_t)NULL;
}

static void handle_click_wait(struct FTR *f, int b, int m, int x, int y)
{
	if (b == 1 || b == 2 || b == 3)
		ftr_notify_the_desire_to_stop_this_loop(f, 10000*y + x);
}

void ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y, int *b, int *m)
{
	ftr_set_handler(f, "button", handle_click_wait);
	int r = ftr_loop_run(f);
	if (x) *x = r % 10000;
	if (y) *y = r / 10000;
}


// vim:set foldmethod=marker:
