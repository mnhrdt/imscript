// FTR simple code examples
//
// The FTR library provides windows with the following properties:
// 	* the pixels can be accessed by a user-editable framebuffer
// 	* key presses, mouse clicks and movement can be captured
//
// Note: FTR stands for "Finestra-Teclat-Ratolí", which is catalan for
// "Window-Keyboard-Mouse".  It is slightly more pronounceable than WKM.
//

#include "ftr.h"


// minimal program: open a black window, and exit the program when it is closed
int main_minimal()
{
	return ftr_loop_run(ftr_new_window());
}

// the same program as above, but without a memory leak
int main_minimal2()
{
	struct FTR *f = ftr_new_window();
	int r = ftr_loop_run(f);
	free(f);
	return r;
}

// minimal (alternative API)
int main_minimal()
{
	return ftr_loop_run(&ftr_new_window());
}

// minimal2 (alternative API)
int main_minimal2()
{
	struct FTR f = ftr_new_window();
	return ftr_loop_run(&f);
}

// minimal3 (new window with drawable buffer)
int main_minimal3()
{
	int w = 320;
	int h = 200;
	struct FTR f = ftr_new_window_and_malloc_buffer(w, h);
}

// minimal4 (image display)
int main_minimal3()
{
	// load image
	int w, h;
	char *buf = load_image_rgb24("file.png", &w, &d);

	// create window
	struct FTR f = ftr_new_window_with_image_rgb24(buf, w, h);

	return ftr_loop_run(&f);
}



// open an image and display it inside a window
int main_viewimage()
{
	// load image
	int w, h;
	char *buf = load_image_rgb("file.png", &w, &d);

	// create window
	struct FTR f = ftr_new_window();

	// associate image to window
	ftr_put_image_in_window(&f, buf, w, h);

	// run
	return ftr_loop_run(&f);
}

// like above, but with fancier features
int main_viewimage2(void)
{
	// create window
	struct FTR *f = ftr_new_window();

	// load image
	f->buf = load_image_rgb("file.png", &f->w, &f->h);
	f->changed = 1;

	// handler to close window when ESC is pressed
	// (this handler is set by default)
	//ftr_set_handler(f, "key", ftr_handler_exit_on_ESC);

	// open an image window on a different process
	ftr_loop_fork(f);

	// do other stuff
	printf("i'm here!\n");

	// when finished, close the other process (optional)
	ftr_quit(f);

	return 0;
}

// like above, but with fancier features
int main_viewimage2_alt(void)
{
	// load image
	int w, h;
	char *buf = load_image_rgb("file.png", &w, &h);

	// create window
	struct FTR f = ftr_new_window();

	// handler to close window when ESC is pressed
	//ftr_set_handler(&f, "key", ftr_handler_exit_on_ESC);

	// open an image window on a different process
	ftr_loop_fork(&f);

	// do other stuff
	printf("i'm here!\n");

	// when finished, close the other process (optional)
	ftr_quit(f);

	return 0;
}

// like above, but using a convenience wrapper
int main_viewimage3(void)
{
	// load image
	int w, h;
	char *buf = load_image_rgb("file.png", &w, &h);

	ftr_spawn_image(f, buf, w, h);

	printf("i'm here!\n");
	free(buf); // note: the image is still visible and unchangeable!

	return 0;
}
// note: after quitting this program, the forked window continues to exist

// a "fire" demo
int main_fire(void)
{
	struct FTR *f = ftr_new_window();
	ftr_set_handler(f, "idle", my_fire_iteration);
	ftr_set_handler(f, "key", ftr_handler_exit_on_ESC);
	return ftr_loop_run(f);
}


// a "crop" program
int main_crop(void)
{
	// load image
	int w, h;
	char *buf = load_image_rgb("-", &w, &h);

	struct FTR f = ftr_new_window_with_image(buf, w, h);

	int from[2], to[2];
	ftr_wait_for_mouse_click(&f, from+0, from+1, NULL);
	ftr_wait_for_mouse_click(&f, to+0, to+1, NULL);

	ftr_close(&f);

	do_crop(buf, &w, &h, from, to);

	save_image_rgb("-", buf, w, h);

	free(buf);
	return 0;
}




// Two modes of interaction:
// 1. procedural (actions explicitly called by the programmer)
// 2. event-based (actions triggered by events: mouse  click, screen refresh)
//


// pensem-ho en abstracte... ¿com ens agradaria que fos?

// single procedural call
int main_spawn_image()
{
	// load image
	int w, h, pd;
	char *buf = load_image_rgb("file.png", &w, &f, &pd);

	// show image in its own window
	ftr_spawn_image(buf, w, h, pd);

	// continue with program
	free(buf);
	return 0;
}

// Requirement: FTR *never* allocates space for image data
// Data always and only belongs to the user

// single procedural call
int main_spawn_image()
{
	// load image
	int w, h;
	char *buf = load_image_rgb24("file.png", &w, &h);

	// show image in its own window
	ftr_new_window_with_image_rgb24(buf, w, h);

	// continue with program
	free(buf);
	return 0;
}

// multiple procedural calls
int main_spawn_image2()
{
	// load image
	int w, h;
	char *buf = load_image_rgb24("file.png", &w, &h);

	// show image in its own window
	struct FTR f = ftr_new_window_with_image_rgb24(buf, w, h);

	// set first pixel to black
	buf[0] = buf[1] = buf[2] = 0;

	// tell window that the buffer has changed
	ftr_update(&f);

	// do other stuff
	sleep(10);

	// close window
	ftr_close(&f);

	// continue with program
	free(buf);
	return 0;
}

// loop of procedural calls
int main_fire2()
{
	int w = 320, h = 200;
	float *buf = malloc(w*h*sizeof*buf);
	memset(buf, 0, w*h*sizeof*buf);

	struct FTR f = ftr_new_window_with_image_float_gray(buf, w, h);

	for (int i = 0 ; i < 2000; i++)
	{
		fire_iteration(buf, w, h);
		ftr_update(&f);
	}

	ftr_close(&f);

	return 0;
}
// observation: the user of this program can close the window before
// the 2000 iterations are done.  Nevertheless, the remaining iterations
// are computed (and the calls to ftr_update fail silently).  The next
// program is slightly more aware.

// loop of procedural calls with exit control
int main_fire2()
{
	int w = 320, h = 200;
	float *buf = malloc(w*h*sizeof*buf);
	memset(buf, 0, w*h*sizeof*buf);

	struct FTR f = ftr_new_window_with_image_float_gray(buf, w, h);

	while(ftr_is_open(&f)) {
		fire_iteration(buf, w, h);
		ftr_update(&f);
	}

	free(buf);
	return 0;
}

// loop of procedural calls with blocking events
int main_polygon()
{
	int w = 320, h = 200;
	char *buf = malloc(3*w*h);
	memset(buf, 0, 3*w*h);

	struct FTR f = ftr_new_window_with_image_rgb24(buf, w, h);

	while(ftr_is_open(&f)) {
		printf("click for new point\n");
		int p, q, b;
		ftr_wait_for_mouse_click(&f, &p, &q, &b);
		printf("got click at (%d %d) with button %d!\n\n", p, q, b);
	}


	return 0;
}


// simplest event-driven program
int main_events()
{
	ftr_loop_run(&ftr_new_window())'
}
