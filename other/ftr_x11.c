//#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
//#include <X11/Xutil.h> // only for XDestroyImage, that can be easily removed
#include <X11/XKBlib.h> // only for XkbKeycodeToKeysym
#include <unistd.h> // only for "fork"
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

	// X11-only internal data
	Display *display;
	Visual *visual;
	Window window;
	GC gc;
	XImage *ximage;
	int imgupdate;
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];


// for debug purposes
char *event_names[] ={
"Nothing		0",
"None			1",
"KeyPress		2",
"KeyRelease		3",
"ButtonPress		4",
"ButtonRelease		5",
"MotionNotify		6",
"EnterNotify		7",
"LeaveNotify		8",
"FocusIn		9",
"FocusOut		10",
"KeymapNotify		11",
"Expose			12",
"GraphicsExpose		13",
"NoExpose		14",
"VisibilityNotify	15",
"CreateNotify		16",
"DestroyNotify		17",
"UnmapNotify		18",
"MapNotify		19",
"MapRequest		20",
"ReparentNotify		21",
"ConfigureNotify	22",
"ConfigureRequest	23",
"GravityNotify		24",
"ResizeRequest		25",
"CirculateNotify	26",
"CirculateRequest	27",
"PropertyNotify		28",
"SelectionClear		29",
"SelectionRequest	30",
"SelectionNotify	31",
"ColormapNotify		32",
"ClientMessage		33",
"MappingNotify		34",
"GenericEvent		35",
"LASTEvent		36"};



struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];

	// general
	f->w = w;
	f->h = h;
	f->max_w = 2000;
	f->max_h = 2000;
	f->rgb = malloc(f->max_w * f->max_h * 3);
	for (int i = 0; i < 3*w*h; i++)
		f->rgb[i] = x ? x[i] : 0;

	// specific to X11 backend
	f->display = XOpenDisplay(NULL);
	if (!f->display) exit(fprintf(stderr, "Cannot open display\n"));
	int s = DefaultScreen(f->display);
	int white = WhitePixel(f->display, s);
	int black = BlackPixel(f->display, s);
	f->gc = DefaultGC(f->display, s);
	f->visual = DefaultVisual(f->display, s);
	f->window = XCreateSimpleWindow(f->display,
			RootWindow(f->display, s), 10, 10, f->w, f->h, 1,
			black, white);
			//white, black);
	f->ximage = NULL;
	f->imgupdate = 1;
	XSelectInput(f->display, f->window, (1L<<25)-1-(1L<<7)-ResizeRedirectMask);
	//XSelectInput(f->display, f->window, ExposureMask | KeyPressMask);
	//XSelectInput(f->display, f->window, (1L<<25)-1);
	//XSelectInput(f->display, f->window, (1L<<25)-1);
	//XSelectInput(f->display, f->window,
	//	       	ExposureMask
	//		| KeyPressMask
	//		| ButtonPressMask
	//		| PointerMotionMask
	//		//| ResizeRedirectMask
	//		| StructureNotifyMask
	//		);
	XMapWindow(f->display, f->window);

	// general again
	f->handle_key = ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->stop_loop = 0;
	f->changed = 0;

	// run the loop until the first expose
	f->handle_expose2 = ftr_handler_stop_loop;
	ftr_loop_run((struct FTR *)f);
	f->handle_expose2 = 0;

	return *(struct FTR *)f;
}

struct FTR ftr_new_window(int w, int h)
{
	return ftr_new_window_with_image_uint8_rgb(NULL, w, h);
}

void ftr_close(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	if (f->ximage) f->ximage->f.destroy_image(f->ximage);
	if (f->rgb) free(f->rgb);
	XCloseDisplay(f->display);
}

static int keycode_to_ascii(struct _FTR *f, int keycode, int keystate)
{
	//fprintf(stderr, "\tkeycode = %d (%d)\n", keycode, keystate);
	//int key = XKeycodeToKeysym(f->display, keycode, keystate);
	int key = XkbKeycodeToKeysym(f->display, keycode, 0, keystate);
	if (keycode == 9)   return 27;    // ascii ESC
	if (keycode == 119) return 127;   // ascii DEL
	if (keycode == 22)  return '\b';
	if (keycode == 23)  return '\t';
	if (keycode == 36 || keycode == 105) return '\n';
	return key;
}

static int do_bound(int a, int b, int x)
{
	if (b < a) return do_bound(b, a, x);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

//static int do_snap(int a, int b, int x)
//{
//	if (b < a) return do_snap(b, a, x);
//	if (2*x > a+b) return b;
//	return a;
//}
//
//static void do_doublesnap(int *ox, int *oy, int w, int h, int x, int y)
//{
//	if (x < 0) x = 0;
//	if (y < 0) y = 0;
//	if (x >= w) x = w-1;
//	if (y >= h) y = h-1;
//	int dhor = fmin(x, w - x);
//	int dver = fmin(y, h - y);
//	if (dhor < 100 && dver > 100) { x = x<100 ? 0 : w-1; }
//	if (dver < 100 && dhor > 100) { y = y<100 ? 0 : h-1; }
//	*ox = x;
//	*oy = y;
//}

static void process_next_event(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	XEvent event;
	XNextEvent(f->display, &event);
	//fprintf(stderr, "ev %d\t\"%s\"\n", event.type,
	//		event_names[event.type]);

	if (event.type == Expose || f->changed) {
		if (f->handle_expose)
			f->handle_expose(ff, 0, 0, 0, 0);
		f->changed = 0;

		if (!f->ximage || f->imgupdate) {
			if (f->ximage) f->ximage->f.destroy_image(f->ximage);
			f->ximage = XGetImage(f->display, f->window,
					0, 0, f->w, f->h, AllPlanes, ZPixmap);
			f->imgupdate = 0;
		}
		for (int i = 0; i < f->w*f->h; i++) {
			f->ximage->data[4*i+0] = f->rgb[3*i+2];
			f->ximage->data[4*i+1] = f->rgb[3*i+1];
			f->ximage->data[4*i+2] = f->rgb[3*i+0];
			f->ximage->data[4*i+3] = -1;
		}
		XPutImage(f->display, f->window, f->gc, f->ximage,
				0, 0, 0, 0, f->w, f->h);
		if (f->handle_expose2)
			f->handle_expose2(ff, 0, 0, 0, 0);

	}

	if (event.type == ConfigureNotify) {
		XConfigureEvent e = event.xconfigure;
		if (f->w != e.width || f->h != e.height) {
			f->w = e.width < f->max_w ? e.width : f->max_w;
			f->h = e.height< f->max_h ? e.height : f->max_h;
			if (f->handle_resize)
				f->handle_resize(ff, 0, 0, f->w, f->h);
			f->imgupdate = 1;
		}
	}


	// call the motion handler also for enter/leave events
	// (this is natural for most applications)
	if (event.type == EnterNotify && f->handle_motion) {
		XCrossingEvent e = event.xcrossing;
		//fprintf(stderr, "enter notify (%d %d)\n", e.x, e.y);
		f->handle_motion(ff, -1, e.state, e.x, e.y);
	}
	if (event.type == LeaveNotify && f->handle_motion) {
		XCrossingEvent e = event.xcrossing;
		int x = do_bound(0, f->w - 1, e.x);
		int y = do_bound(0, f->h - 1, e.y);
		//fprintf(stderr,"LEAVE notify (%d %d)[%d %d]\n",e.x,e.y,x,y);
		f->handle_motion(ff, -1, e.state, x, y);
	}

	// "expose" and "resize" are never lost
	// the mouse and keyboard evends are ignored when too busy
	if (XPending(f->display) > 0) return;


	if (event.type == ButtonPress && f->handle_button) {
		XButtonEvent e = event.xbutton;
		f->handle_button(ff, e.button, e.state, e.x, e.y);
	}
	if (event.type == KeyPress && f->handle_key) {
		XKeyEvent e = event.xkey;
		//f->handle_key(ff, e.keycode, e.state, e.x, e.y);
		//int keysym = XKeycodeToKeysym(f->display, e.keycode, e.state);
		int key = keycode_to_ascii(f, e.keycode, e.state);
		f->handle_key(ff, key, e.state, e.x, e.y);
	}
	if (event.type == MotionNotify && f->handle_motion) {
		XMotionEvent e = event.xmotion;
		f->handle_motion(ff, e.is_hint, e.state, e.x, e.y);
	}
}

int ftr_loop_run(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	while (!f->stop_loop)
	{
		if (!f->handle_idle || XPending(f->display) > 0)
			process_next_event(ff);
		else if (f->handle_idle) {
			XEvent ev;
			ev.type = Expose;
			XLockDisplay(f->display);
			XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
			XFlush(f->display);
			XUnlockDisplay(f->display);
			f->handle_idle(ff, 0, 0, 0, 0);
			f->changed = 1; // not sure if necessary...
		}
	}

	int r = f->stop_loop;
	f->stop_loop = 0;
	return r;
}

int ftr_num_pending(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	return XPending(f->display);
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
		f->stop_loop = 1;
}

void ftr_handler_exit_on_ESC_or_q(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || tolower(k)=='q')
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
	if (0) ;
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
	if (0) ;
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
