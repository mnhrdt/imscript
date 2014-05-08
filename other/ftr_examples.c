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

// open an image and display it inside a window
int main_viewimage()
{
	// load image
	int w, h;
	char *buf = load_image_rgb("file.png", &w, &f);

	// create window
	struct FTR *f = ftr_new_window();

	// associate image to window
	f->w = w;
	f->h = h;
	f->buf = buf;
	f->changed = 1;

	// run
	return ftr_loop_run(f);
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
	ftr_set_handler(f, "key", ftr_handler_exit_on_ESC);

	// open an image window on a different process
	ftr_loop_fork(f);

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
