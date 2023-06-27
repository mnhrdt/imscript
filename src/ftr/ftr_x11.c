//#define _POSIX_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <signal.h>
#include <X11/Xlib.h>
//#include <X11/Xutil.h> // only for XDestroyImage, that can be easily removed
#include <unistd.h> // only for "fork"

#include <assert.h>


#include "ftr.h"

// X11-specific part {{{1
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
	ftr_event_handler_t handle_wheel;
	int max_w, max_h;
	int stop_loop;

	// X11-only internal data
	Display *display;
	Visual *visual;
	Window window;
	GC gc;
	XImage *ximage;
	int imgupdate;

	int wheel_ax;

	pid_t child_pid;

	int s;    // display scaling, read from envvar at startup
	int W, H; // scaled sizes, for convenience
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


static int getenv_int(char *s, int d)
{
	char *t = getenv(s);
	return t ? atoi(t) : d;
}

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

	// scaling hack
	f->s = getenv_int("FTR_SCALING", 1);
	f->W = f->w * f->s;
	f->H = f->h * f->s;
	//fprintf(stderr, "w,h = %d,%d\n", f->w, f->h);
	//fprintf(stderr, "s = %d\n", f->s);
	//fprintf(stderr, "W,H = %d,%d\n", f->W, f->H);

	// specific to X11 backend
	f->display = XOpenDisplay(NULL);
	if (!f->display) exit(fprintf(stderr, "Cannot open display\n"));
	int s = DefaultScreen(f->display);
	int white = WhitePixel(f->display, s);
	int black = BlackPixel(f->display, s);
	f->gc = DefaultGC(f->display, s);
	f->visual = DefaultVisual(f->display, s);
	f->window = XCreateSimpleWindow(f->display, RootWindow(f->display, s),
			10, 10, f->W, f->H, 1, black, white);
			//white, black);
	XStoreName(f->display, f->window, "ftr");
	f->ximage = NULL;
	f->imgupdate = 1;
	f->wheel_ax = 0;
	int mask = 0
		| ExposureMask
		| KeyPressMask
		| ButtonPressMask
		| ButtonReleaseMask
		| PointerMotionMask
		| StructureNotifyMask
		| EnterWindowMask
		| LeaveWindowMask
		;
	XSelectInput(f->display, f->window, mask);
	XMapWindow(f->display, f->window);

	// immanentize the pid
	{
		Atom a = XInternAtom(f->display, "_NET_WM_PID", 0);
		int p = getpid();
		XChangeProperty(f->display, f->window,
				a, 6, 32, 0, (void*)&p, 1);
	}

	// general again
	f->handle_key = ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->handle_wheel = NULL;
	f->stop_loop = 0;
	f->changed = 0;

	// run the loop until the first expose
	f->handle_expose2 = ftr_handler_stop_loop;
	ftr_loop_run((struct FTR *)f);
	f->handle_expose2 = 0;

	//return *(struct FTR *)f;
	// avoid stupid type-punning warning:
	struct FTR F[1];
	memcpy(F, f, sizeof f);
	return *F;
}

void ftr_change_title(struct FTR *ff, char *s)
{
	struct _FTR *f = (void*)ff;
	XStoreName(f->display, f->window, s);
}

void ftr_close(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	if (f->ximage) f->ximage->f.destroy_image(f->ximage);
	if (f->rgb) free(f->rgb);
	XCloseDisplay(f->display);
}

// replacement for XKeycodeToKeysym, which is deprecated
static int x_keycode_to_keysym(struct _FTR *f, int keycode)
{
	int nothing;
	KeySym *t = XGetKeyboardMapping(f->display, keycode, 1, &nothing);
	int r = t[0];
	XFree(t);
	return r;
}

static int keycode_to_ftr(struct _FTR *f, int keycode, int keystate)
{
	int key = x_keycode_to_keysym(f, keycode);
	//fprintf(stderr, "keycode to keysym : %d => %d\n", keycode, key);

	if (keycode == 9)   return 27;    // ascii ESC
	if (keycode == 119) return 127;   // ascii DEL
	if (keycode == 22)  return '\b';
	if (keycode == 23)  return '\t';
	if (keycode == 36 || keycode == 105) return '\n';
	if (keycode == 111) return FTR_KEY_UP;
	if (keycode == 113) return FTR_KEY_LEFT;
	if (keycode == 114) return FTR_KEY_RIGHT;
	if (keycode == 116) return FTR_KEY_DOWN;
	if (keycode == 112) return FTR_KEY_PAGE_UP;
	if (keycode == 117) return FTR_KEY_PAGE_DOWN;
	return key;
}

static int do_bound(int a, int b, int x)
{
	if (b < a) return do_bound(b, a, x);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static void ugly_hack_disable_enter_events(struct _FTR *f)
{
	int mask = 0
		| ExposureMask
		| KeyPressMask
		| ButtonPressMask
		| ButtonReleaseMask
		| PointerMotionMask
		| StructureNotifyMask
		;
	XSelectInput(f->display, f->window, mask);
}

static void process_next_event(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	XEvent event = { .type = GenericEvent } ;
	if (!f->changed)
		XNextEvent(f->display, &event); // blocks and waits!
	//fprintf(stderr,"ev(%p,%p) %d\t\"%s\"\n",(void*)f->display,(void*)f->window,event.type,event_names[event.type]);

	if (event.type == Expose || f->changed) {
		if (f->handle_expose)
			f->handle_expose(ff, 0, 0, 0, 0);
		f->changed = 0;

		if (!f->ximage || f->imgupdate) {
			if (f->ximage) f->ximage->f.destroy_image(f->ximage);
			//f->ximage = XShmGetImage(f->display, f->window,
			//		0, 0, f->w, f->h, AllPlanes, ZPixmap);
			f->ximage = XGetImage(f->display, f->window,
				0, 0, f->W, f->H, AllPlanes, ZPixmap);
			f->imgupdate = 0;
		}
		if (f->s == 1)
		for (int i = 0; i < f->w*f->h; i++) {
			// drain bramage
			f->ximage->data[4*i+0] = f->rgb[3*i+2];
			f->ximage->data[4*i+1] = f->rgb[3*i+1];
			f->ximage->data[4*i+2] = f->rgb[3*i+0];
			f->ximage->data[4*i+3] = 0;
		}

		// OK, boys, this is where the actual fucking happens
		if (f->s > 1) {
			int s = f->s;
			for (int j = 0; j < f->h; j++)
			for (int i = 0; i < f->w; i++)
			{
				int ij = i + j*f->w;
				int IJ = s*i + s*j*f->W;
				for (int q = 0; q < s; q++)
				for (int p = 0; p < s; p++)
				{
					int Ip = IJ + p + q*f->W;
					f->ximage->data[4*Ip+0]= f->rgb[3*ij+2];
					f->ximage->data[4*Ip+1]= f->rgb[3*ij+1];
					f->ximage->data[4*Ip+2]= f->rgb[3*ij+0];
					f->ximage->data[4*Ip+3]= 0;
				}
			}
		}
		//XShmPutImage(f->display, f->window, f->gc, f->ximage,
		//		0, 0, 0, 0, f->w, f->h, 0);
		XPutImage(f->display, f->window, f->gc, f->ximage,
				0, 0, 0, 0, f->W, f->H);
		if (f->handle_expose2)
			f->handle_expose2(ff, 0, 0, 0, 0);

	}

	if (event.type == ConfigureNotify) {
		XConfigureEvent e = event.xconfigure;
		if (f->W != e.width || f->H != e.height) {
			f->w= e.width/f->s  < f->max_w ? e.width/f->s:f->max_w;
			f->h= e.height/f->s < f->max_h ? e.height/f->s:f->max_h;
			f->W = f->w * f->s;
			f->H = f->h * f->s;
			if (f->handle_resize)
				f->handle_resize(ff, 0, 0, f->w, f->h);
			f->imgupdate = 1;
		}
	}

	if (event.type == ButtonPress && f->handle_wheel) {
		XButtonEvent e = event.xbutton;
		if (e.button == 4 || e.button == 5) {
			int np = XPending(f->display);
			f->wheel_ax += e.button==4 ? 1 : -1;
			fprintf(stderr, "\twheel acc(%d) %d\n", f->wheel_ax,np);
			if (np < 2) {
				int x = e.x / f->s;
				int y = e.y / f->s;
				f->handle_wheel(ff, f->wheel_ax, e.state, x,y);
				f->wheel_ax = 0;
			}
			return;
		}
	}

	// "expose" and "resize" are never lost
	// the mouse and keyboard events are ignored when too busy
	//
	int ne = XPending(f->display);
	if (ne > 1) // magical "1" here
	//if (ne > 0) // magical "1" here
	//{
	//	fprintf(stderr, "\tlost{%d} %s\n", ne, event_names[event.type]);
	      return;
	//}

	// call the motion handler also for enter/leave events
	// (this is natural for most applications)
	if (!f->handle_wheel) {//wheel accumulation does not work well with this
	if (event.type == EnterNotify && f->handle_motion) {
		XCrossingEvent e = event.xcrossing;
		//fprintf(stderr, "enter notify (%d %d)\n", e.x, e.y);
		if (XPending(f->display)) return;
		f->handle_motion(ff, -1, e.state, e.x / f->s, e.y / f->s);
	}
	if (event.type == LeaveNotify && f->handle_motion) {
		XCrossingEvent e = event.xcrossing;
		int x = do_bound(0, f->w - 1, e.x / f->s);
		int y = do_bound(0, f->h - 1, e.y / f->s);
		//fprintf(stderr,"LEAVE notify (%d %d)[%d %d]\n",e.x,e.y,x,y);
		if (XPending(f->display)) return;
		f->handle_motion(ff, -2, e.state, x, y);
	}
	}


	if (event.type == MotionNotify && f->handle_motion) {
		XMotionEvent e = event.xmotion;
		f->handle_motion(ff, e.is_hint, e.state, e.x/f->s, e.y/f->s);
	}
	if (event.type == ButtonPress && (f->handle_button||f->handle_wheel)) {
		XButtonEvent e = event.xbutton;
		if (f->handle_wheel && (e.button == 4 || e.button == 5)) {
			f->wheel_ax += e.button==4 ? 1 : -1;
			fprintf(stderr, "\twheel ack %d\n", f->wheel_ax);
			//if (!XPending(f->display)) {
				int x = e.x / f->s;
				int y = e.y / f->s;
				f->handle_wheel(ff, f->wheel_ax, e.state, x,y);
				f->wheel_ax = 0;
			//}
		}
		f->handle_button(ff, 1<<(7+e.button),e.state,e.x/f->s,e.y/f->s);
		static int ugly_disable = 1;
		if (ugly_disable && (e.button == 4 || e.button == 5))
		{
			ugly_hack_disable_enter_events(f);
			ugly_disable = 0;
		}
	}
	if (event.type == ButtonRelease && f->handle_button) {
		XButtonEvent e = event.xbutton;
		if (e.button != 4 && e.button != 5)
			f->handle_button(ff, -(1<<(7+e.button)),
					e.state, e.x/f->s, e.y/f->s);
	}
	if (event.type == KeyPress && f->handle_key) {
		XKeyEvent e = event.xkey;
		//f->handle_key(ff, e.keycode, e.state, e.x, e.y);
		//int keysym = XKeycodeToKeysym(f->display, e.keycode, e.state);
		if (e.keycode != 77) {
		int key = keycode_to_ftr(f, e.keycode, e.state);
		f->handle_key(ff, key, e.state, e.x/f->s, e.y/f->s);
		}
	}
}

void ftr_x11_force_redraw(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	XEvent e;
	e.type = Expose;
	XSendEvent(f->display, f->window, 0, NoEventMask, &e);
	XFlush(f->display);
}

int ftr_loop_run(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;

	//f->changed = 1;

	while (!f->stop_loop)
	{
		if (!f->handle_idle || f->changed || XPending(f->display) > 0)
			process_next_event(ff);
		else if (f->handle_idle) {
			f->handle_idle(ff, 0, 0, 0, 0);
			XEvent ev;
			ev.type = Expose;
			//XLockDisplay(f->display);
			XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
			XFlush(f->display);
			//XUnlockDisplay(f->display);
		}
	}

	int r = f->stop_loop;
	f->stop_loop = 0;
	return r;
}

int ftr_loop_run2(struct FTR *ff, struct FTR *gg)
{
	struct _FTR *f = (void*)ff;
	struct _FTR *g = (void*)gg;

	char *dn_f = XDisplayString(f->display);
	char *dn_g = XDisplayString(g->display);
	if (0 != strcmp(dn_f, dn_g))
		exit(fprintf(stderr, "FTR error: two displays bad bad bad (%p,%p)(\"%s\",\"%s\")\n", (void*)f->display, (void*)g->display, dn_f, dn_g));

	fprintf(stderr, "dn_f = %p \"%s\"\n", (void*)f->display, dn_f);
	fprintf(stderr, "dn_g = %p \"%s\"\n", (void*)g->display, dn_g);

	while (!f->stop_loop && !g->stop_loop)
	{
		int pending_f = XPending(f->display);
		int pending_g = XPending(g->display);

		fprintf(stderr, "pending fg = %d %d\n", pending_f, pending_g);

		// treat g
		if (!g->handle_idle || g->changed || pending_g > 0)
			process_next_event(gg);
		else if (g->handle_idle) {
			g->handle_idle(gg, 0, 0, 0, 0);
			XEvent ev;
			ev.type = Expose;
			XSendEvent(g->display, g->window, 0, NoEventMask, &ev);
			XFlush(g->display);
		}

		// treat f
		if (!f->handle_idle || f->changed || pending_f > 0)
			process_next_event(ff);
		else if (f->handle_idle) {
			f->handle_idle(ff, 0, 0, 0, 0);
			XEvent ev;
			ev.type = Expose;
			XSendEvent(f->display, f->window, 0, NoEventMask, &ev);
			XFlush(f->display);
		}
	}

	int r = f->stop_loop + g->stop_loop;
	f->stop_loop = g->stop_loop = 0;
	return r;
}

int ftr_num_pending(struct FTR *ff)
{
	struct _FTR *f = (void*)ff;
	return XPending(f->display);
}

void ftr_notify_the_desire_to_stop_this_loop(struct FTR *ff, int retval)
{
	struct _FTR *f = (void*)ff;
	f->stop_loop = retval?retval:-1;
}

// common to all implementations {{{1


#include "ftr_common_inc.c"

// vim:set foldmethod=marker:
