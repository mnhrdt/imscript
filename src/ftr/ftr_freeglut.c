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
	int w, h;
	unsigned char *rgb;
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
	int max_w, max_h;
	int stop_loop;

	// glut-only internal data
	// (there's none, glut requires using global variables)
	//int glut_stopped;
	int handle_mute;
	int glut_initted;
	int glut_button_mask;
	int glut_keymod_mask;
	int glut_window_x;
	int glut_window_y;

	// pixel scaling (internal)
	int s;
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];

void ftr_change_title(struct FTR *ff, char *s)
{
}

// glut-specific part {{{1

// global variable
struct _FTR *ftr_freeglut_global_state = NULL;

static void my_idle(void);

// glut handlers {{{2
static void my_displayfunc(void)
{
	//fprintf(stderr, "GLUT displayfunc\n");

	struct _FTR *f = ftr_freeglut_global_state;

	if (f->handle_expose && f->changed) {
		f->handle_expose((void*)f, 0, 0, 0, 0);
		f->changed = 0;
	}

	if (f->handle_idle)
		glutIdleFunc(my_idle); // ugly hack to allow toggling

	glViewport(0, 0, f->w, f->h);

	glRasterPos2i(-1,1);
	glPixelZoom(f->s, -f->s);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // ugly

	//fprintf(stderr, "glDrawPixels %p\n", f->rgb);
	glDrawPixels(f->w, f->h, GL_RGB, GL_UNSIGNED_BYTE, f->rgb);

	glutSwapBuffers();
}

static void my_reshapefunc(int w, int h)
{
	//fprintf(stderr, "GLUT reshapefunc %d %d\n", w, h);

	struct _FTR *f = ftr_freeglut_global_state;

	f->w = w < f->max_w ? w : f->max_w;
	f->h = h < f->max_h ? h : f->max_h;

	if (f->handle_resize)
		f->handle_resize((void*)f, 0, 0, f->w, f->h);
}

static int key_from_glut_to_ftr(int k)
{
	return k + 1000;
}

static int button_from_glut_to_ftr(int b)
{
	return 1 << ( b + 8 );
}

static int modifiers_from_glut_to_ftr(int m)
{
	int r = 0;
	if (m & GLUT_ACTIVE_SHIFT) r |= FTR_MASK_SHIFT;
	if (m & GLUT_ACTIVE_CTRL ) r |= FTR_MASK_CONTROL;
	if (m & GLUT_ACTIVE_ALT  ) r |= FTR_MASK_MOD1;
	return r;
}

static void my_specialfunc(int k, int x, int y)
{
	//fprintf(stderr, "GLUT specialfunc  %d '%c' %d %d\n",
	//		k, isalnum(k)?k:' ', x, y);

	struct _FTR *f = ftr_freeglut_global_state;
	if (f->handle_key) {
		int m = modifiers_from_glut_to_ftr(glutGetModifiers());
		f->handle_key((void*)f, key_from_glut_to_ftr(k), m,
				x/f->s, y/f->s);
		if (f->changed)
			glutPostRedisplay();
	}
}

static void my_keyboardfunc(unsigned char k, int x, int y)
{
	//fprintf(stderr, "GLUT keyboardfunc  %d '%c' %d %d\n",
	//		k, isalnum(k)?k:' ', x, y);

	struct _FTR *f = ftr_freeglut_global_state;

	if (f->handle_key) {
		int m = modifiers_from_glut_to_ftr(glutGetModifiers());
		f->handle_key((void*)f, k, m, x/f->s, y/f->s);
		if (f->changed)
			glutPostRedisplay();
	}
}

static void my_mousefunc(int b, int s, int x, int y)
{
	//fprintf(stderr, "GLUT mousefunc %d %d (%d %d)\n", b, s, x, y);

	struct _FTR *f = ftr_freeglut_global_state;

	int B = button_from_glut_to_ftr(b);
	if (s == 0) // pressed
		f->glut_button_mask |= B;
	if (s == 1) // de-pressed
		f->glut_button_mask &= ~B;


	if (f->handle_button /*&& 0 == s*/) {
		if (s && (B == FTR_BUTTON_UP || B == FTR_BUTTON_DOWN))
			return;
		int m = modifiers_from_glut_to_ftr(glutGetModifiers());
		f->handle_button((void*)f, B*(s?-1:1), m, x/f->s, y/f->s);
		if (f->changed)
			glutPostRedisplay();
	}
}

static void my_passivemotionfunc(int x, int y)
{
	struct _FTR *f = ftr_freeglut_global_state;
	if (f->handle_mute) return;

	//fprintf(stderr, "GLUT passive motionfunc %d %d\n", x, y);

	f->handle_mute = 1;
	if (f->handle_motion) {
		int m = modifiers_from_glut_to_ftr(glutGetModifiers());
		if (f->changed)
			glutPostRedisplay();
		f->handle_motion((void*)f, 0, m|f->glut_button_mask,
				x/f->s, y/f->s);
	}
	f->handle_mute = 0;
}

static void my_motionfunc(int x, int y)
{
	//fprintf(stderr, "GLUT motionfunc %d %d\n", x, y);
	struct _FTR *f = ftr_freeglut_global_state;
	my_passivemotionfunc(x, y);
}


static void my_idle(void)
{
	//fprintf(stderr, "GLUT idle handler\n");

	struct _FTR *f = ftr_freeglut_global_state;
	if (f->handle_idle) {
		f->handle_idle((void*)f, 0, 0, 0, 0);
		if (f->changed)
			glutPostRedisplay();
	} else {
		glutIdleFunc(NULL);
	}
}

static char *name_with_pid(void)
{
	static char n[FILENAME_MAX];
	pid_t p = getpid();
	snprintf(n, FILENAME_MAX, "ftr_win_pid_%d", p);
	return n;
}

// setup_glut_environment {{{2
static void setup_glut_environment(struct _FTR *f)
{
	int argc = 2;
	char *argv[]={"dummy", "-gldebug", NULL};
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
			GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutCreateWindow(name_with_pid());
	glutReshapeWindow(f->w * f->s, f->h * f->s);
	glutDisplayFunc(my_displayfunc);
	glutReshapeFunc(my_reshapefunc);
	glutIdleFunc(NULL);
	glutMouseFunc(my_mousefunc);
	glutMotionFunc(my_motionfunc);
	glutPassiveMotionFunc(my_passivemotionfunc);
	glutKeyboardFunc(my_keyboardfunc);
	glutSpecialFunc(my_specialfunc);
	f->glut_initted = 1;
	f->glut_button_mask = 0; // not sure
	f->glut_keymod_mask = 0;
	//fprintf(stderr, "setup glut rgb = %p\n", f->rgb);
	if (f->glut_window_x >= 0)
		glutPositionWindow(f->glut_window_x, f->glut_window_y);
}

static int getenv_int(char *s, int d)
{
	char *t = getenv(s);
	return t ? atoi(t) : d;
}


// ftr_new_window_with_image_uint8_rgb {{{2
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];

	f->w = w;
	f->h = h;
	f->max_w = 4000;
	f->max_h = 3000;
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
	f->glut_button_mask = 0;
	f->glut_window_x = f->glut_window_y = -1;
	f->s = getenv_int("FTR_SCALING", 1);

	// freeglut specific part
	setup_glut_environment(f);

	//f->handle_expose2 = ftr_handler_stop_loop;
	//ftr_loop_run((struct FTR *)f);
	//f->handle_expose2 = 0;

	//return *(struct FTR *)f;
	// avoid stupid type-punning warning:
	struct FTR F[1];
	memcpy(F, f, sizeof f);
	return *F;
}

// ftr_close {{{2
void ftr_close(struct FTR *ff)
{
	// do nothing by now

	//struct _FTR *f = (void*)ff;
	//if (f->ximage) XDestroyImage(f->ximage);
	//if (f->rgb) free(f->rgb);
	//XCloseDisplay(f->display);
}

// ftr_loop_run {{{2
int ftr_loop_run(struct FTR *ff)
{
	//fprintf(stderr, "going to start loop\n");
	struct _FTR *f = (void*)ff;

	// after window creation, we jumped into the air
	// here, we get hold of a new branch
	ftr_freeglut_global_state = f;

	if (!f->glut_initted)
		setup_glut_environment(f);

	if (f->handle_idle)
		glutIdleFunc(my_idle);

	glutMainLoop();
	//fprintf(stderr, "returned from glut loop (%d)\n", f->stop_loop);

	int r = f->stop_loop;
	f->stop_loop = 0;
	return r;
}

// ftr_notify_the_desire_to_stop_this_loop {{{2
void ftr_notify_the_desire_to_stop_this_loop(struct FTR *ff, int retval)
{
	//fprintf(stderr, "stop notification retval = %d\n", retval);
	struct _FTR *f = (void*)ff;
	f->stop_loop = retval;
	f->glut_initted = 0;

	// save window position
	f->glut_window_x = glutGet(GLUT_WINDOW_X);
	f->glut_window_y = glutGet(GLUT_WINDOW_Y);

	glutLeaveMainLoop();
}

int ftr_x11_force_redraw(struct FTR *ff)
{
	// just ignore it
	(void)ff;
}

#include "ftr_common_inc.c"

// vim:set foldmethod=marker:
