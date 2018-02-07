---
layout: netc
title: FTR, a minimal window-keyboard-mouse library
---
<h1>FTR</h1>

<p>
FTR is a portable and ported library for opening a <b>window</b>, drawing
pixels to it, and reading <b>keyboard</b> and <b>mouse</b> input.  It is
crippled, slow, and extremely minimal.  The main goal is to write tiny
portable programs that open windows and show images.  It is a thin wrapper
around native libraries, intended to be as thin as possible, but not thinner.
</p>

<p>
Current ports include X11 (tested on linux, freebsd and openbsd),
ANSI-terminal (without mouse),  and any platform supporting <a
href="http://freeglut.sourceforge.net/">freeglut</a> (tested on linux, should
work on Windows and OSX after installing freeglut).  Planned ports include
native Windows in the short term and native OSX on the long term.
</p>

<h2>Examples</h2>

<h3>mini.c</h3>

<p>
This is the simplest program using FTR:
{% highlight c %}
#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	return ftr_loop_run(&f);
}
{% endhighlight %}
</p>

<p>
To compile and run this program on a unix system with X, place the file "ftr.c"
besides "mini.c" and run
{% highlight sh %}
cc mini.c -lX11
./a.out
{% endhighlight %}
</p>

<p>
This should open a black window of size 320&times;200 that can be closed by
pressing ESC.
</p>

<h3>events.c</h3>

<p>
A useful program is obtained by capturing some events on the window:
{% highlight c %}

#include <stdio.h>
#include "ftr.c"

void print_event(struct FTR *f, int k, int m, int x, int y)
{
	printf("event k=%d m=%d x=%d y=%d\n", k, m, x, y);
}

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	ftr_set_handler(&f, "key", print_event);
	ftr_set_handler(&f, "button", print_event);
	ftr_set_handler(&f, "motion", print_event);
	ftr_set_handler(&f, "resize", print_event);
	return ftr_loop_run(&f);
}
{% endhighlight %}
</p>

<p>
This program opens a black window and prints a line when an "event" happens
on that window.  Events include key-presses, mouse motions and clicks, and
window resizing.  The key-presses are represented by the ASCII value of the
corresponding key (using the correct keyboard configuration).
</p>

<h3>greenpixel.c</h3>

<p>
The following program opens a black window with a green pixel at the position
(10,10):
{% highlight c %}
#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	f.rgb[1 + 3 * (f.w * 10 + 10)] = 255;
	f.changed = 1;
	return ftr_loop_run(&f);
}
{% endhighlight %}
</p>

<p>
First, notice that the <tt>struct FTR</tt> is not an opaque data type, but it
has data that can be read and written explicitly, such as the framebuffer
<tt>rgb</tt> and the image size.
</p>

<p>
Second, observe the line <tt>f.changed = 1</tt>, that tells "the system" that
the framebuffer contents have changed.  This line is necessary
because, as opposed to many other "widget" toolkits, the call to
<tt>ftr_new_window</tt> actually opens the window.  The call to
<tt>ftr_loop_run</tt> is optional, and only needed so that the program does
not exit immediately.  To understand that, try comenting the line
<tt>f.changed=1</tt>.  The window will be opened but the green pixel will not
appear immediately; it will appear only when the window is redrawn, such as
when moving another window in front of this one.  Another way to understand
this is to run the following program:
{% highlight c %}
#include <unistd.h>
#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	return sleep(3);
}
{% endhighlight %}
</p>

<p>
This program opens a window for three seconds, and then it exists and closes
the window.  This observation is very important: FTR is not an event-based
library, it is multi-paradigm.  You can write perfectly useful programs
without a main loop, or with many main loops in different places.  There are
API functions like <tt>ftr_wait_for_mouse_click</tt> that allow a synchronous
programming style, which is impossible in other toolkits such as GLUT:
{% highlight c %}
#include <stdio.h>
#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	int x, y;
	ftr_wait_for_mouse_click(&f, &x, &y);
	printf("the mouse was clicked at position %d %d\n", x, y);
	return 0;
}
{% endhighlight %}
</p>

<h3>fire.c</h3>

<p>
The "idle" event is called continuously.  This allows, for example, to run
simulations in realtime.  However, it is dangerous because the CPU is used
100%.  The following program draws a fire effect (compile with optimizations
to obtain an acceptable framerate)
{% highlight c %}
#include <stdio.h>
#include <stdlib.h>
#include "ftr.c"

static void draw_fire(struct FTR *f, int x, int y, int k, int m)
{
	int num_lines_bottom = 3;

	// build palette
	unsigned char palette[3*256];
	for (int i = 0; i < 256; i++) {
		palette[3 * i + 0] = 4 * i;
		palette[3 * i + 1] = (255-2 * i)/3;
		palette[3 * i + 2] = 3 * i;
	}

	// build buffer (if resize, restart buffer)
	static float *t = NULL;
	static int w = 0;
	static int h = 0;
	if (!f || w != f->w || h != f->h) { 
		w = f->w;
		h = f->h;
		if (t) free(t);
		t = malloc(w * h * sizeof*t); 
		for (int i = 0; i < w*h; i++)
			t[i] = (unsigned char)rand();
	}

	// draw random values at the bottom
	int p = 0;
	for (int j = 0; j < num_lines_bottom; j++)
	for (int i = 0; i < w; i++) {
		t[p] = (unsigned char)(t[p] + 15*(rand()/(1.0+RAND_MAX)));
		p++;
	}

	// paint pixels by combining lower rows
	for (int j = h-1; j >= num_lines_bottom; j--)
	for (int i = 0; i < w; i++) {
		p = j*w+i;
		t[p] = (t[p-w] + 2 * t[p-2*w-1] + 2 * t[p-2*w] + 2 * t[p-2*w+1])
			* 9 * 4 / 256;
	}

	// render with palette
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = (unsigned char)(t[w*(h-j-1) + i]);
		f->rgb[3 * (w*j+i) + 0] = palette[3 * idx + 0];
		f->rgb[3 * (w*j+i) + 1] = palette[3 * idx + 1];
		f->rgb[3 * (w*j+i) + 2] = palette[3 * idx + 2];
	}
}

int main(void)
{
	struct FTR f = ftr_new_window(800, 600);
	ftr_set_handler(&f, "idle", draw_fire);
	return ftr_loop_run(&f);
}
{% endhighlight %}

<canvas id="myCanvas" width="320" height="200">caca</canvas>

<script type="text/javascript">
var c = document.getElementById("myCanvas");
var ctx = c.getContext("2d");
var id = ctx.getImageData(0, 0, c.width, c.height);
var w = id.width;
var h = id.height;

// Build palette.
var palette = Array();
for (var i = 0; i < 256; i++) {
  palette[3 * i + 0] = (4 * i) % 256;
  palette[3 * i + 1] = (2 * i) % 256;
  palette[3 * i + 2] = (1 * i) % 256;
}

var f = Array();
f.length = w * h;
for (var i = 0; i < f.length; i++) f[i] = Math.random() * 255;

function draw() {
  var num_lines_bottom = 3;

  // Draw random values at the bottom.
  var p = w * (h - num_lines_bottom);
  for (var i = 0; i < num_lines_bottom; i++)
    for (var j = 0; j < w; j++) {
      f[p] = (f[p] + 15 * Math.random()) % 256;
      p++;
    }

  // Paint pixels by combining lower rows.
  p = 0;
  for (var i = 0; i < h - num_lines_bottom; i++) {
    for (var j = 0; j < w; j++) {
      f[p] = (f[p+w] + 2 * f[p+2*w-1] + 2 * f[p+2*w] + 2 * f[p+2*w+1])
           * 9 * 4 / 256;
      p++;
    }
  }

  // Render with palette.
  p = 0;
  for (var i = 0; i < h; i++)
    for (var j = 0; j < w; j++) {
      var idx = Math.floor(f[p]);
      id.data[4 * p + 0] = palette[3 * idx + 0];
      id.data[4 * p + 1] = palette[3 * idx + 1];
      id.data[4 * p + 2] = palette[3 * idx + 2];
      id.data[4 * p + 3] = 255;
      p++;
    }

  ctx.putImageData(id, 0, 0);
}

setInterval("draw()", 50);

</script>
(Javascript fire by Pau Gargallo)
</p>

<h3>ftr.h</h3>

<p>
Besides the examples above, the header file <tt>ftr.h</tt> is a useful and
complete source of information:
{% highlight c %}
#ifndef _FTR_H
#define _FTR_H

// data structure to store the state of a window
struct FTR {
	// visible state
	int w, h;           // size of the image
	unsigned char *rgb; // rgb 24-bit image
	int changed;        // variable to indicate that the image has changed

	void *userdata;     // ignored by the library

	// hidden implementation details
	char pad[100];
};

// type of a handler function
typedef void (*ftr_event_handler_t)(struct FTR*,int,int,int,int);

// core API
struct FTR ftr_new_window(int w, int h);
void       ftr_change_title(struct FTR *f, char *title);
int        ftr_set_handler(struct FTR *f, char *event, ftr_event_handler_t h);
int        ftr_loop_run(struct FTR *f);
void       ftr_notify_the_desire_to_stop_this_loop(struct FTR *f, int retval);

// convenience functions (reducible to the core api)
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *i, int w, int h);
void       ftr_close(struct FTR *f);
void       ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y);
void       ftr_wait_for_key_depress(struct FTR *f, int *key, int *modifiers);
int        ftr_num_pending(struct FTR *f);
void       ftr_loop_fork(struct FTR *f);
void       ftr_fork_window_with_image_uint8_rgb(unsigned char *i, int w, int h);
ftr_event_handler_t ftr_get_handler(struct FTR *f, char *id);

// example handlers
void ftr_handler_exit_on_ESC     (struct FTR*,int,int,int,int);
void ftr_handler_exit_on_ESC_or_q(struct FTR*,int,int,int,int);
void ftr_handler_toggle_idle     (struct FTR*,int,int,int,int);
void ftr_handler_stop_loop       (struct FTR*,int,int,int,int);
void ftr_handler_dummy           (struct FTR*,int,int,int,int);


// ascii keys with name (necessary because '\e' is not standard)
#define FTR_KEY_ESC        27
#define FTR_KEY_DEL        127

// non-ascii keys (numbers inspired by glut)
#define FTR_KEY_FN         1000
#define FTR_KEY_LEFT       1100
#define FTR_KEY_UP         1101
#define FTR_KEY_RIGHT      1102
#define FTR_KEY_DOWN       1103
#define FTR_KEY_PAGE_UP    1104
#define FTR_KEY_PAGE_DOWN  1105
#define FTR_KEY_HOME       1106
#define FTR_KEY_END        1107
#define FTR_KEY_INSERT     1108

// key modifiers (numbers inspired by X)
#define FTR_MASK_SHIFT     1
#define FTR_MASK_LOCK      2
#define FTR_MASK_CONTROL   4
#define FTR_MASK_MOD1      8
#define FTR_MASK_MOD2      16
#define FTR_MASK_MOD3      32
#define FTR_MASK_MOD4      64
#define FTR_MASK_MOD5      128

// buttons and button modifiers (inspired by X)
#define FTR_BUTTON_LEFT    256
#define FTR_BUTTON_MIDDLE  512
#define FTR_BUTTON_RIGHT   1024
#define FTR_BUTTON_UP      2048
#define FTR_BUTTON_DOWN    4096

#endif//_FTR_H
{% endhighlight %}
</p>

<h2>Principles</h2>

<h3>Goals of FTR</h3>

<ul>
	<li>Be written in standard C</li>
	<li>No global state in the API</li>
	<li>The "ftr.h" file is as small as possible (ideally, one printed
	page)</li>
	<li>The "ftr.c" file contains the whole library</li>
	<li>The "ftr.c" file can be directly included, so the header is
	unnecessary</li>
	<li>Deal only with WINDOWS, KEYBOARD and MOUSE</li>
	<li>The keyboard input is read in ascii, when possible</li>
	<li>The mouse wheel works, somehow</li>
	<li>It works "natively" in unix, windows and mac</li>
	<li>The same API hides different backends: Xlib, windows, cocoa,
	glut, text</li>
	<li>In Unix, it uses only Xlib (already implemented)</li>
	<li>In Windows, it uses the Windows Api (not yet implemented)</li>
	<li>In mac, it uses the COCOA C interface (not yet implemented)</li>
	<li>There is a fall-back GLUT backend (not yet implemented)</li>
	<li>There is a fall-back freeglut backend (already implemented)</li>
	<li>There is a fall-back text-based backend, where events are read
	from stdin</li>
	<li>Programs must have exactly the same behaviour on all
	backends</li>
	<li>The complete API is implemented on all supported backends</li>
	<li>Multiple user interface paradigms are supported:
	<ol>
		<li>event-based, automatic loop (set event handlers and call
		ftr_loop)</li>
		<li>event based, explicit loop (write the event loop
		yourself)</li>
		<li>blocking calls (explicitly wait for user input)</li>
		<li>forked windows (each window has its own process, and
		communication is done via signals and shared memory)</li>
		<li>contiuously-redrawn window showing the evolution of an
		image</li>
	</li>
</ol>
<li>A pixel is defined as three bytes containing the RGB components.</li>
<li>An image is an array of pixels</li>
<li>Do not use typedef structs, useless casts or otherwise ugly C code</li>
<li>80 columns, tabs are 8 spaces, functions have less than 25 lines</li>
<li>To port FTR to a new unix backend, you have to write only 7 "core"
functions:<br />
	<tt>ftr_new_window</tt><br />
	<tt>ftr_change_title</tt><br />
	<tt>ftr_set_handler</tt><br />
	<tt>ftr_loop_run</tt><br />
	<tt>ftr_notify_that_the_framebuffer_has_changed</tt><br />
	<tt>ftr_notify_the_desire_to_stop_this_loop</tt><br />
	<tt>ftr_close_window</tt><br />
	</li>
<li>To port FTR to a non-unix backend, this other function has to be
re-written:<br />
	<tt>ftr_fork_window</tt><br />
	</li>
</ul>

<h3>Non-goals of FTR</h3>

<ul>
<li>Be written in other languages than C</li>
<li>Image input/output</li>
<li>Any kind of image processing</li>
<li>Drawing widgets (e.g., buttons)</li>
<li>Drawing text</li>
<li>Drawing primitives (segments, disks, ...)</li>
<li>OpenGL access</li>
<li>Accelerated or otherwise efficient operations</li>
<li>Deal with anything else than WINDOWS, KEYBOARD or MOUSE</li>
<li>Pixel types other than 24bit RGB</li>
<li>Transparency, broken pixels (separate color channels), gray-level images,
etc</li>
<li>Floating-point samples</li>
</ul>

<p>
The API may change to adapt to future needs, and new functions may be added
without problem.  However, <b>the non-goals are strong decisions and I accept no
	compromise on them</b>.  For example, FTR will never perform a contrast change.
This is up to the user of the library.  On the other hand, some "helper" and
"convenience" functions are defined to ease common usages.  For example, the
"ftr_open_window_with_image" function.
</p>

<h3>Status of ports</h3>

<p>The X11 interface is tested and working in linux and mac.
The freeglut interface should be working in any platform that has glut, but
needs some testing, specially in the keyboard handler.</p>

<table border="1" cellpadding="6" cellspacing="0">
	<tr>
		<th>Back-end</th>
		<th>Type</th>
		<th>Status</th>
	</tr>
	<tr>
		<td>X11</th>
		<td>native</td>
		<td>Working, tested in linux and OSX--Quartz</td>
	</tr>
	<tr>
		<td>cocoa</th>
		<td>native</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>windows</th>
		<td>native</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>glut</th>
		<td>semi-native</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>glfw</th>
		<td>semi-native</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>freeglut</th>
		<td>semi-native</td>
		<td>Working, tested in linux (some keyboard issues)</td>
	</tr>
	<tr>
		<td>gdk</th>
		<td>experimental</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>fltk</th>
		<td>experimental</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>qt</th>
		<td>experimental</td>
		<td><i>not written</i></td>
	</tr>
	<tr>
		<td>xcb</th>
		<td>experimental</td>
		<td><i>not written</i></td>
	</tr>
</table>

<h3>Motivation</h3>

<p>
FTR is inspired and tries to imitate many other libraries.  The principal
sources of inspiration are <a href="http://cimg.sourceforge.net/">CImg</a>
and <a href="http://www.opengl.org/resources/libraries/glut/">GLUT</a>.
</p>


<h4>Why not use CImg?</h4>

<p>
CImg is essentially perfect.
There's no real reason to not use CImg in C++.  Sometimes I would
like to use CImgDdisplay only, without need for CImg.  For example, when I want
to publish an image processing algorithm of my own (that does not use the CImg
image processing capabilities), and I want to provide a small interface.  CImg
is perfect for that, but a bit overkill if you do not work with the CImg class.
CImg.h is 46.000 lines of code.  My objective is that ftr.c should never go
beyond 1000.
</p>


<h4>Why not use GLUT?</h4>

<p>
GLUT is essentially perfect, it has a very beautiful API
and provides a neat cross-platform way to open a window and show pixels on it.
There are, however a few problems: you can not exit easily the main loop, the
mouse wheel does not generate events, and you are forced to draw your pixels
using OpenGL.  The freeglut re-implementation of GLUT solves the problem of the
mouse wheel (but is not native on mac osx) and it lets you exit the event loop
(but it closes the window when you do that).  Freeglut is about 20.000 lines of
code in about 30 files, including font data and widget-drawing stuff.
</p>


<h4>How is FTR different?</h4>

<p>
Besides its limited scope (no widgets, no image processing, no opengl
access), FTR gives you more flexibility in the way you write your program.
The idea is that you write your whole program (e.g., an algorithm in image
processing), without needing FTR at all.  Then, as an afterthought, you
decide to add, at some scattered places in your program, a few calls to FTR
to do things like:
</p>

<p> 1. Nested deep inside your algorithm, you want to open a window
	containing a temporary image, and then continue running the  program.
	The window is kept open until it is closed by the user or until the
	program ends.  This involves a single call to this function:<br />
	<br />

	<tt>ftr_new_window_with_image(char *rgb, int w, int h);</tt><br />
	<br />
	and nothing else.  No initialization, no boilerplate.  This function
	returns immediately.</p>

<p>2. In the main loop of your program, you want to see how the iterations
progress.  Thus, you open a window with the image before the loop, and
inside the loop you update the pixels of the image and then call the
function <tt>ftr_notify_that_the_framebuffer_has_changed</tt>.</p>

<p>3. Your example program uses a hard-coded pixel position.  By calling
<tt>ftr_wait_for_mouse_click</tt> you let the user to click on a pixel.</p>

<p> 4. See the tutorial for further examples.</p>
</ol>


<h2>Download, etc.</h2>

<p>
FTR is under heavy development and has not yet released a distribution.
By "distribution" I mean a single file <tt>ftr.c</tt>.  Right now, I have two
files
<tt><a href="https://github.com/mnhrdt/imscript/blob/master/other/ftr_x11.c">ftr_x11.c</a></tt>
<tt><a
		href="https://github.com/mnhrdt/imscript/blob/master/other/ftr_freeglut.c">ftr_freeglut.c</a></tt>,
containing the respective implementations, that help me understand the required
symmetries.
</p>

<p>
You can find these files inside the "other" subdirectory of <a
	href="https://github.com/mnhrdt/imscript/">imscript</a>
</p>


<h2>Prize</h2>

<p>
The first user who guesses correcly the meaning of the acronym "FTR" will be
invited to a beer of his choice.
</p>

<div style="text-align:right;font-size:small;color:gray;">
	<i>Last updated: 4 March 2016, Enric Meinhardt-Llopis</i>
</div>


<!-- vim:set tw=77 filetype=html spell spelllang=en: -->
