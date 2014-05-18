// c99 -O3 xhw.c -lX11

#define _POSIX_C_SOURCE 199309L
#include <assert.h>
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
	assert(t->tv_sec >= first_seconds);
	double r = (t->tv_sec - first_seconds) + 1e-9*t->tv_nsec;
	return r;
}

#include <X11/Xlib.h>
#include <X11/Xutil.h> // only for XDestroyImage, that can be easily removed
#include <X11/XKBlib.h> // only for XkbKeycodeToKeysym

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <ctype.h>


void firechar(unsigned char *buf, int w, int h);


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

	int changed;

	// X11-related internal data
	Display *display;
	Visual *visual;
	Window window;
	GC gc;
	XImage *ximage;
	int imgupdate;
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
	//f->w = 1000;
	//f->h = 1000;
	f->max_w = 1000;
	f->max_h = 1000;
	int bufsize = 4 * f->max_w * f->max_h;
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
	f->imgupdate = 1;

	f->changed = 0;

	//XSelectInput(f->display, f->window, ExposureMask | KeyPressMask);
	XSelectInput(f->display, f->window, (1L<<25)-1-(1L<<7)-ResizeRedirectMask);
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

	f->handle_key = NULL;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->do_exit = 0;

	return f;
}

int ftr_loop_run(struct FTR *f)
{
	char *msg = "Hello, World!";
	XEvent event;


	while (!f->do_exit)
	{
		if (XPending(f->display) > 0) {
		XNextEvent(f->display, &event);
		//fprintf(stderr, "ev %d\t\"%s\"\n", event.type,
		//		event_names[event.type]);

		//static int pos = 0;
		//for (int i = 0; i < 42; i++)
		//	f->rgb8_buffer[pos++] = 0;

		if (event.type == Expose || f->changed)
		{
			f->changed = 0;
			//fprintf(stderr, "\texpose event\n");
			//XFillRectangle(f->display, f->window, f->gc,
			//		20, 20, f->w - 40, f->h - 40);
			//XDrawString(f->display, f->window, f->gc,
			//		50, 12, msg, strlen(msg));

			if (!f->ximage && event.type == Expose) {
				f->ximage = XGetImage(f->display, f->window,
					0, 0, f->w, f->h, AllPlanes, ZPixmap);
			}
			if (f->ximage) {
				if (f->imgupdate) {
				f->ximage->data = f->rgb8_buffer;
				f->ximage->width = f->w;
				f->ximage->height = f->h;
				f->ximage->bytes_per_line = 0;
				if (!XInitImage(f->ximage))
					exit(fprintf(stderr,"e:xinit image\n"));
				f->imgupdate = 0;
				}
				XPutImage(f->display, f->window, f->gc, f->ximage, 0, 0, 0, 0, f->w, f->h);
			}

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
				f->imgupdate = 1;
			}
		}
	}//Xpending
		//fprintf(stderr, "idle\n");
		XLockDisplay(f->display);
		{
			event.type = Expose;
			XSendEvent(f->display, f->window, 0,NoEventMask, &event);
			XFlush(f->display);
		}
		XUnlockDisplay(f->display);
	firechar((unsigned char *)f->rgb8_buffer, f->w, f->h);
	f->changed = 1;
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

void firechar(unsigned char *buf, int w, int h)
{
	static int cx = 0;
	cx += 1;
	if (0 == cx % 100) {
		double static os = 0;
		double s = seconds();
		double dif = s - os;
		double fps = 100/dif;
		fprintf(stderr, "CX = %d (%g) {%g} %f FPS\n", cx, s, dif, fps);
		os = s;
	}
	static unsigned char *pal = NULL;
	if (!pal) {
		// build palette
		pal = malloc(3*256);
		for (int i = 0; i < 256; i++) {
			pal[3*i+0] = (1*i)%256;     //blue
			pal[3*i+1] = (3*i)%256;     //blue
			pal[3*i+2] = (4*i)%256;     //blue
			//pal[3*i+0] = (2*i)%256;     //blue
			//pal[3*i+1] = i<128?0:2*(i-128); //green
			//pal[3*i+2] = i<128?2*i:255;//i<128?2*i:255;     //red
			//pal[3*i+1] = i<64?0:2*(i-64);
			//pal[3*i+2] = i<64?4*i:255;
		}
		//for (int i = 128; i <= 255; i++)
		//for (int l = 0; l < 3; l++)
		//	pal[3*i+l] = pal[3*(255-i)+l];
	}

	static float *f = NULL;
	static int iw = 0;
	static int ih = 0;
	if (!f || iw!=w || ih!=h) { iw = w; ih = h; f = malloc(w*h*sizeof*f); }
	f[0] = 0;

	//if (cx > 2000) return;
	int num_lines_bottom = 3;

	// draw random values at the bottom
	int p = 0;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		//f[p] = sin(0.1*(f[p] + 0.1*(cos(3.14*(i-w/2.0)/w))*sin(rand())));
		//f[p] += 0.3*(cos(3.14*(i-w/2.0)/w))*sin(rand());
		f[p] = fmod(f[p] + 15*(rand()/(1.0+RAND_MAX)),256);
		p++;
	}

	// paint pixels by combining lower rows
	//for (int j = num_lines_bottom; j < h; j++)
	for (int j = h-1; j >= num_lines_bottom-1; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		f[p] = (1.5*f[p-w] + 1.7 * f[p-2*w+1] + 1.5 * f[p-2*w] + 2.3 * f[p-2*w-1])
			/7.01;
	}

	// render with palette
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int iidx = w*(h-j-1) + i;
		int idx = (unsigned char)(f[iidx]/1.2);
		int pos = w*j + i;
		buf[4*pos+0] = pal[3*idx+0];
		buf[4*pos+1] = pal[3*idx+1];
		buf[4*pos+2] = pal[3*idx+2];
		buf[4*pos+3] = 255;
	}

	////for (int i = 0; i < w*h*3; i++)
	////	buf[i] = rand();
	//int numlines = 3;
	//for (int j = h - 1; j >= h-numlines-1; j--)
	//for (int i = 1; i < w; i++)
	//for (int l = 0; l < 4; l++)
	//	buf[(j*w+i)*4+l] = (buf[(j*w+i-1)*4+l]+(13.0*(rand()/(1.0+RAND_MAX))));

	//for (int j = h-numlines-1; j >= 0; j--)
	//for (int i = 0; i < w; i++)
	//for (int l = 0; l < 4; l++)
	//{
	//	float x1 = buf[((j+1)*w+i-1)*4+l];
	//	float x2 = buf[((j+1)*w+i  )*4+l];
	//	float x3 = buf[((j+1)*w+i+1)*4+l];
	//	float g = 0.7*x2 + 0.2*x1 + 0.1*x3;
	//	buf[(j*w+i)*4+l] = (buf[(j*w+i)*4+l] + g)/2;
	//}
}


// begin the actual program here

void my_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "GOT key %d, %d at (%d %d)\n", k, m, x, y);
	//int ks = XKeycodeToKeysym(f->display, k, 0);
	int ks = XkbKeycodeToKeysym(f->display, k, 0, 0);
	int ks1 = XkbKeycodeToKeysym(f->display, k, 0, m);
	if (isprint(ks))
		fprintf(stderr, "\t'%c' '%c'\n", ks, ks1);
	if (k == 9)
		f->do_exit = 1;
	//fprintf(stderr, "\tw h = %d %d\n", f->w, f->h);
	firechar((unsigned char *)f->rgb8_buffer, f->w, f->h);
	f->changed = 1;
}

void my_button_handler(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "GOT button %d, %d at (%d %d)\n", b, m, x, y);
	firechar((unsigned char *)f->rgb8_buffer, f->w, f->h);
	f->changed = 1;
}

void my_motion_handler(struct FTR *f, int h, int m, int x, int y)
{
	fprintf(stderr, "GOT motion %d, %d at (%d %d)\n", h, m, x, y);
	firechar((unsigned char *)f->rgb8_buffer, f->w, f->h);
	f->changed = 1;
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
	return ftr_loop_run(f);
}




#if 0
void my_draw_frame(FTR *f)
{
	// fes una iteració del foc
}

int main_foc(void)
{
	struct FTR *f = ftr_new_window(320, 200);
	ftr_set_handler(f, "idle", my_draw_frame);
	ftr_set_handler(f, "key", ftr_key_exit_on_ESC);
	int r = ftr_start_loop(f);
	return r;
}



struct FTR *f = my_spawn_image(float *x, int w, int h, int pd)
{
	struct FTR *f = ftr_new_window(w, h);
	ftr_fork_loop(f); // runs on a parallel thread, waiting
	// to be exited or stopped from inside
}

int main_show_image(void)
{
	int w, h, pd;
	float *x = iio_read_image_float_vec("-", &w, &h, &pd);
	struct FTR *f = my_spawn_image(x, w, h, pd);
	...
	...
	...
	ftr_stop(f);
	return 0;
}
#endif






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
