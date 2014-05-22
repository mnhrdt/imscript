
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

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
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];



void my_display(void)
{
	glRasterPos2i(-1,1);
	glPixelZoom(1,-1);

	int w = global_state->pw;
	int h = global_state->ph;
	//unsigned char *x = global_state->view;
	//glDrawPixels(w, h, GL_RGB, GL_UNSIGNED_BYTE, x);
	glDrawPixels(w, h, GL_RGB, GL_FLOAT, global_state->view);

	glutSwapBuffers();
}

struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];

	f->w = w;
	f->h = h;
	f->max_w = 2000;
	f->max_h = 2000;
	f->rgb = malloc(f->max_w * f->max_h * 3);
	for (int i = 0; i < 3*w*h; i++)
		f->rgb[i] = x[i];
	f->changed = 0;

	// freeglut specific part
	int argc = 1;
	char *argv[]={"dummy"};
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutCreateWindow("Whatever");
	glutReshapeWindow(f->w, f->h);
	glutDisplayFunc(my_display);
	glutReshapeFunc(NULL);
	glutIdleFunc(NULL);
	glutMouseFunc(NULL);
	glutMotionFunc(NULL);
	glutPassiveMotionFunc(NULL);
	glutKeyboardFunc(NULL);
	glutSpecialFunc(NULL);

	// X11-specific part
	//f->display = XOpenDisplay(NULL);
	//if (!f->display) exit(fprintf(stderr, "Cannot open display\n"));
	//int s = DefaultScreen(f->display);
	//int white = WhitePixel(f->display, s);
	//int black = BlackPixel(f->display, s);
	//f->gc = DefaultGC(f->display, s);
	//f->visual = DefaultVisual(f->display, s);
	//f->window = XCreateSimpleWindow(f->display, RootWindow(f->display, s),
	//		10, 10, f->w, f->h, 1, black, white);
	//f->ximage = NULL;
	//f->imgupdate = 1;
	//XSelectInput(f->display,f->window,(1L<<25)-1-(1L<<7)-ResizeRedirectMask);
	//XMapWindow(f->display, f->window);




	f->handle_key = ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->stop_loop = 0;

	f->handle_expose2 = ftr_handler_stop_loop;
	ftr_loop_run((struct FTR *)f);
	f->handle_expose2 = 0;

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

//static void process_next_event(struct FTR *ff)
//{
//	struct _FTR *f = (void*)ff;
//
//	glutMainLoopEvent();
//	//XEvent event;
//	//XNextEvent(f->display, &event);
//	////fprintf(stderr, "ev %d\t\"%s\"\n", event.type,
//	////		event_names[event.type]);
//
//	//if (event.type == Expose || f->changed) {
//	//	f->changed = 0;
//
//	//	if (!f->ximage || f->imgupdate) {
//	//		if (f->ximage) XDestroyImage(f->ximage);
//	//		f->ximage = XGetImage(f->display, f->window,
//	//				0, 0, f->w, f->h, AllPlanes, ZPixmap);
//	//		f->imgupdate = 0;
//	//	}
//	//	for (int i = 0; i < f->w*f->h; i++) {
//	//		f->ximage->data[4*i+0] = f->rgb[3*i+2];
//	//		f->ximage->data[4*i+1] = f->rgb[3*i+1];
//	//		f->ximage->data[4*i+2] = f->rgb[3*i+0];
//	//		f->ximage->data[4*i+3] = -1;
//	//	}
//	//	XPutImage(f->display, f->window, f->gc, f->ximage,
//	//			0, 0, 0, 0, f->w, f->h);
//	//	if (f->handle_expose2)
//	//		f->handle_expose2(ff, 0, 0, 0, 0);
//
//	//}
//	//if (event.type == KeyPress && f->handle_key) {
//	//	XKeyEvent e = event.xkey;
//	//	f->handle_key(ff, e.keycode, e.state, e.x, e.y);
//	//}
//	//if (event.type == ButtonPress && f->handle_button) {
//	//	XButtonEvent e = event.xbutton;
//	//	f->handle_button(ff, e.button, e.state, e.x, e.y);
//	//}
//	//if (event.type == MotionNotify && f->handle_motion) {
//	//	XMotionEvent e = event.xmotion;
//	//	f->handle_motion(ff, e.is_hint, e.state, e.x, e.y);
//	//}
//	//if (event.type == ConfigureNotify && f->handle_resize) {
//	//	XConfigureEvent e = event.xconfigure;
//	//	if (f->w != e.width || f->h != e.height) {
//	//		f->w = e.width < f->max_w ? e.width : f->max_w;
//	//		f->h = e.height< f->max_h ? e.height : f->max_h;
//	//		f->handle_resize(ff, 0, 0, f->w, f->h);
//	//		f->imgupdate = 1;
//	//	}
//	//}
//}

int ftr_loop_run(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	while (!f->stop_loop)
		glutMainLoopEvent();
//	{
//
//		if (!f->handle_idle || XPending(f->display) > 0)
//			process_next_event(ff);
//		else if (f->handle_idle) {
//			XEvent ev;
//			ev.type = Expose;
//			XLockDisplay(f->display);
//			XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
//			XFlush(f->display);
//			XUnlockDisplay(f->display);
//			f->handle_idle(ff, 0, 0, 0, 0);
//			f->changed = 1;
//		}
//	}

	int r = f->stop_loop;
	f->stop_loop = 0;
	return r;
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
	if  (k == 9)
		f->stop_loop = 1;
}

void ftr_handler_stop_loop(struct FTR *f, int k, int m, int x, int y)
{
	f->stop_loop = 1;
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
	if (false) ;
	else if (0 == strcmp(id, "key"))    { f->handle_key    = e; return 0; }
	else if (0 == strcmp(id, "button")) { f->handle_button = e; return 0; }
	else if (0 == strcmp(id, "motion")) { f->handle_motion = e; return 0; }
	else if (0 == strcmp(id, "expose")) { f->handle_expose = e; return 0; }
	else if (0 == strcmp(id, "resize")) { f->handle_resize = e; return 0; }
	else if (0 == strcmp(id, "idle"))   { f->handle_idle   = e; return 0; }
	else return fprintf(stderr, "WARNING: unrecognized event \"%s\"\n", id);
}

ftr_event_handler_t ftr_get_handler(struct FTR *ff, char *id)
{
	struct _FTR *f = (void*)ff;
	if (false) ;
	else if (0 == strcmp(id, "key"))    return f->handle_key   ;
	else if (0 == strcmp(id, "button")) return f->handle_button;
	else if (0 == strcmp(id, "motion")) return f->handle_motion;
	else if (0 == strcmp(id, "expose")) return f->handle_expose;
	else if (0 == strcmp(id, "resize")) return f->handle_resize;
	else if (0 == strcmp(id, "idle"))   return f->handle_idle  ;
	else return fprintf(stderr, "WARNING: bad event \"%s\"\n", id),
		(ftr_event_handler_t)NULL;
}

static void handle_click_wait(struct FTR *f, int b, int m, int x, int y)
{
	if (b == 1 || b == 2 || b == 3)
		f->stop_loop = 10000*y + x;
}

void ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y, int *b, int *m)
{
	ftr_set_handler(f, "button", handle_click_wait);
	int r = ftr_loop_run(f);
	if (x) *x = r % 10000;
	if (y) *y = r / 10000;
}
