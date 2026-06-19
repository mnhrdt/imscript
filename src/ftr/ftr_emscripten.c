// NOTE: the emscripten backend is a dry hack, just intended to have some
// pixels showing on a browser and nothing more.  It cannot exit, it cannot
// resize the window, or behave like a civilized browser citizen at all.


#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <emscripten.h>
#include <emscripten/html5.h>

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

	// emscripten specific
	uint8_t *rgba;
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];

// global variable
// (since the emscripten loop does not take any arguments, we need to put
// all state globally accessible)
// 
// this variable is filled-in by the function "ftr_loop_run"
struct _FTR *ftr_emscripten_global_state = NULL;

void em_render(void)
{
	struct _FTR *f = ftr_emscripten_global_state;

	if (!f->changed) return;
	f->changed = 0;

	EM_ASM({
		const c = Module.canvas.getContext('2d');
		const p = $0; // rgba
		const w = $1; // w
		const h = $2; // h
		const v = new Uint8ClampedArray(Module.HEAPU8.buffer, p, w*h*4);
		const i = new ImageData(v, w, h);
		c.putImageData(i, 0, 0);
		},
	f->rgba, f->w, f->h);
}




void ftr_change_title(struct FTR *ff, char *s)
{
	(void)ff;
	(void)s;
}

static void ftr_dump_rgb_to_rgba(struct FTR *ff)
{
	assert(global_w == ff->w);
	assert(global_h == ff->h);
	for (int i = 0; i < f->w * f->h; i++)
	{
		global_rgba[4*i+0] = ff->rgb[3*i+0];
		global_rgba[4*i+1] = ff->rgb[3*i+1];
		global_rgba[4*i+2] = ff->rgb[3*i+2];
		global_rgba[4*i+3] = 255;
	}
	global_dirty = 1;
}


// ftr_new_window_with_image_uint8_rgb {{{2
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];

	f->w = f->max_w = w;
	f->h = f->max_h = h;
	f->rgb = malloc(w * h * 3);
	for (int i = 0; i < 3*w*h; i++)
		f->rgb[i] = x ? x[i] : 0;

	global_w = w;
	global_h = h;
	global_rgba = malloc(w * h * 4);
	global_dirty = 1;
	ftr_dump_rgb_to_rgba(f);

	f->handle_key = NULL;//ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->stop_loop = 0;

	// emscripten shit
	emscripten_set_canvas_element_size("#canvas", global_w, global_h);
	emscripten_set_mousedown_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_mousemove_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_mouseup_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_wheel_callback( "#canvas", NULL, 1, wheel_cb);
	const char *target = EMSCRIPTEN_EVENT_TARGET_DOCUMENT;
	emscripten_set_keydown_callback(target, NULL, 1, key_cb);
	emscripten_set_keyup_callback(target, NULL, 1, key_cb);

	//return *(struct FTR *)f;
	// avoid stupid type-punning warning:
	struct FTR F[1];
	memcpy(F, f, sizeof f);
	return *F;
}

// emscripten handlers
static void my_displayfunc(void)

// ftr_close {{{2
void ftr_close(struct FTR *ff)
{
	fprintf(stderr, "FTR CLOSE!\n");
}

// ftr_loop_run {{{2
int ftr_loop_run(struct FTR *ff)
{
	//fprintf(stderr, "going to start loop\n");
	struct _FTR *f = (void*)ff;

	fprintf(stderr, "entering loop run!\n");
	if (f->handle_expose)
		f->handle_expose(ff, 0, 0, 0, 0);
	ftr_dump_rgb_to_rgba(f);


	// the loop
	emscripten_set_main_loop(render, 0, 1);

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
}

void ftr_x11_force_redraw(struct FTR *ff)
{
	// just ignore it
	(void)ff;
}

#define FTR_NO_FORK
#include "ftr_common_inc.c"

// vim:set foldmethod=marker:
