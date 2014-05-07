/*
 *   * Simple Xlib application drawing a box in a window.
 *     * gcc input.c -o output -lX11
 *       */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>



struct FTR;
typedef void (*ftr_event_handler_t)(struct FTR*,int,int,int,int);


struct FTR {
	// state
	int w, h, max_w, max_h;
	char *rgb8_buffer;
	int do_exit;

	// user-supplied handlers
	ftr_event_handler_t handle_key;
	ftr_event_handler_t handle_button;
	ftr_event_handler_t handle_motion;
	ftr_event_handler_t handle_expose;
	ftr_event_handler_t handle_resize;
	//ftr_event_handler handle_key;
	//ftr_event_handler handle_keydepress;

	// X11-related internal data
	Display *display;
	Visual *visual;
	Window window;
	GC gc;
	XImage *ximage;
};


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

struct FTR *ftr_new_window(void)
{
	struct FTR *f = malloc(sizeof*f);

	f->display = XOpenDisplay(NULL);
	if (!f->display)
		exit(fprintf(stderr, "Cannot open display\n"));

	f->w = 320;
	f->h = 200;
	f->max_w = 1000;
	f->max_h = 1000;
	int bufsize = 3 * f->max_w * f->max_h;
	f->rgb8_buffer = malloc(bufsize);
	memset(f->rgb8_buffer, 0, bufsize);
	for (int i = 0; i < bufsize; i++)
		f->rgb8_buffer[i] = rand();

	int s = DefaultScreen(f->display);
	int white = WhitePixel(f->display, s);
	int black = BlackPixel(f->display, s);
	f->gc = DefaultGC(f->display, s);
	f->visual = DefaultVisual(f->display, s);
	f->window = XCreateSimpleWindow(f->display,
			RootWindow(f->display, s), 10, 10, f->w, f->h, 1,
			black, white);
	f->ximage = NULL;

	//XSelectInput(f->display, f->window, ExposureMask | KeyPressMask);
	//XSelectInput(f->display, f->window, (1L<<25)-1-(1L<<7));
	//XSelectInput(f->display, f->window, (1L<<25)-1);
	//XSelectInput(f->display, f->window, (1L<<25)-1);
	XSelectInput(f->display, f->window,
		       	ExposureMask
			| KeyPressMask
			| ButtonPressMask
			| PointerMotionMask
			//| ResizeRedirectMask
			| StructureNotifyMask
			);

	XMapWindow(f->display, f->window);

	f->handle_key = NULL;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->do_exit = 0;

	return f;
}

int ftr_start_loop(struct FTR *f)
{
	char *msg = "Hello, World!";
	XEvent event;

	while (!f->do_exit)
	{
		XNextEvent(f->display, &event);
		//fprintf(stderr, "ev %d\t\"%s\"\n", event.type,
		//		event_names[event.type]);

		if (event.type == Expose)
		{
			//fprintf(stderr, "\texpose event\n");
			XFillRectangle(f->display, f->window, f->gc,
					20, 20, f->w - 40, f->h - 40);
			XDrawString(f->display, f->window, f->gc,
					50, 12, msg, strlen(msg));

			if (!f->ximage) {
				f->ximage = XGetImage(f->display, f->window,
					0, 0, f->w, f->h, AllPlanes, ZPixmap);
				f->ximage->data = f->rgb8_buffer;
			}
			f->ximage->width = f->w;
			f->ximage->height = f->h;
			f->ximage->bytes_per_line = 0;
			if (!XInitImage(f->ximage))
				exit(fprintf(stderr,"e:xinit image\n"));

			//int s = DefaultScreen(f->display);
			//int white = WhitePixel(f->display, s);
			//int black = BlackPixel(f->display, s);
			//XPutPixel(f->ximage, 1, 1, white);
			//XPutPixel(f->ximage, 2, 2, black);
			XPutImage(f->display, f->window, f->gc, f->ximage, 0, 0, 0, 0, f->w, f->h);
		}
		if (event.type == KeyPress && f->handle_key)
		{
			XKeyEvent e = event.xkey;
			f->handle_key(f, e.keycode, e.state, e.x, e.y);
		}
		if (event.type == ButtonPress && f->handle_button)
		{
			XButtonEvent e = event.xbutton;
			f->handle_button(f, e.button, e.state, e.x, e.y);
		}
		if (event.type == MotionNotify && f->handle_motion)
		{
			XMotionEvent e = event.xmotion;
			f->handle_motion(f, e.is_hint, e.state, e.x, e.y);
		}
		//if (event.type == ResizeRequest && f->handle_resize)
		//{
		//	fprintf(stderr, "WARNING! ResizeRequest\n");
		//	XResizeRequestEvent e = event.xresizerequest;
		//	int rw = e.width  < f->max_w ? e.width  : f->max_w;
		//	int rh = e.height < f->max_h ? e.height : f->max_h;
		//	f->handle_resize(f, 0, 0, rw, rh);
		//}
		if (event.type == ConfigureNotify && f->handle_resize)
		{
			XConfigureEvent e = event.xconfigure;
			if (f->w != e.width || f->h != e.height)
			{
				f->w = e.width < f->max_w ? e.width : f->max_w;
				f->h = e.height< f->max_h ? e.height : f->max_h;
				f->handle_resize(f, 0, 0, f->w, f->h);
			}
		}
	}
	XCloseDisplay(f->display);// actually, not needed
	if (f->ximage) XDestroyImage(f->ximage);
	//free(f->rgb8_buffer);
	free(f);
	return 0;

}

void ftr_framebuffer_dump(struct FTR *f)
{
}

int ftr_set_handler(struct FTR *f, char *id, ftr_event_handler_t e)
{
	if (false) ;
	else if (0 == strcmp(id, "key"))    { f->handle_key    = e; return 0; }
	else if (0 == strcmp(id, "button")) { f->handle_button = e; return 0; }
	else if (0 == strcmp(id, "motion")) { f->handle_motion = e; return 0; }
	else if (0 == strcmp(id, "expose")) { f->handle_expose = e; return 0; }
	else if (0 == strcmp(id, "resize")) { f->handle_resize = e; return 0; }
	else return fprintf(stderr, "WARNING: unrecognized event \"%s\"\n", id);
}



// begin the actual program here

void my_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "GOT key %d, %d at (%d %d)\n", k, m, x, y);
	if (k == 9)
		f->do_exit = 1;
}

void my_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "GOT button %d, %d at (%d %d)\n", b, m, x, y);
}

void my_motion_handler(struct FTR *f, int h, int m, int x, int y)
{
	fprintf(stderr, "GOT motion %d, %d at (%d %d)\n", h, m, x, y);
}

void my_resize_handler(struct FTR *f, int x, int y, int w, int d)
{
	fprintf(stderr, "GOT resize %d %d   %d %d\n", x, y, w, d);
}

int main(void)
{
	struct FTR *f = ftr_new_window();
	ftr_set_handler(f, "key", my_key_handler);
	ftr_set_handler(f, "button", my_button_handler);
	ftr_set_handler(f, "motion", my_motion_handler);
	ftr_set_handler(f, "resize", my_resize_handler);
	return ftr_start_loop(f);
}










// The client creates a connection with the server by calling XOpenDisplay. It
// then requests the creation of a window with XCreateSimpleWindow. A separate
// call to XMapWindow is necessary for mapping the window, that is, for making
// it visible on the screen.
//
// The square is drawn by calling XFillRectangle. This operation can only be
// performed after the window is created. However, performing it once may not
// be enough. Indeed, the content of the window is not always guaranteed to be
// preserved. For example, if the window is covered and then uncovered again,
// its content might require being redrawn. The program is informed that the
// window or a part of it has to be drawn by the reception of an Expose event.
//
// The drawing of the window content is therefore made inside the loop handling
// the events. Before entering this loop, the events the application is
// interested in are selected, in this case with XSelectInput.  The event loop
// waits for an incoming event: if this event is a key press, the application
// exits; if it is an expose event, the window content is drawn. The function
// XNextEvent blocks and flushes the request buffer if there is no event in the
// queue.
