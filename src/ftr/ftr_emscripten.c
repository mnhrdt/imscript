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

	// pixel scaling (internal)
	int s;

	// emscripten specific
	uint8_t *rgba;
	int last_m;
	int last_x;
	int last_y;
};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];

// global variable
// (since the emscripten loop does not take any arguments, we need to put
// all state globally accessible)
// 
// this variable is filled-in before the run, by the function "ftr_loop_run"
struct _FTR *ftr_emscripten_global_state = NULL;


// main emscripten render function
void em_render(void)
{
	struct _FTR *f = ftr_emscripten_global_state;
	if (!f->changed) return;

	if (f->handle_expose)
	{
		//printf("recompute expose\n");
		f->handle_expose((void*)f, 0, 0, 0, 0);
	}

	for (int i = 0; i < f->w * f->h; i++)
	{
		f->rgba[4*i+0] = f->rgb[3*i+0];
		f->rgba[4*i+1] = f->rgb[3*i+1];
		f->rgba[4*i+2] = f->rgb[3*i+2];
		f->rgba[4*i+3] = 255;
	}
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

// documentation for the emscripten/html5 event interface :
// https://emscripten.org/docs/api_reference/html5.h.html
//
// source code (good for macro names) :
// https://github.com/emscripten-core/emscripten/blob/main/system/include/emscripten/html5.h

static int key_from_emkey(const char *s)
{
	if (strlen(s) == 1)
		return *s;
	if (0 == strcmp(s, "PageUp"    )) return FTR_KEY_PAGE_UP;
	if (0 == strcmp(s, "PageDown"  )) return FTR_KEY_PAGE_DOWN;
	if (0 == strcmp(s, "Escape"    )) return FTR_KEY_ESC;
	if (0 == strcmp(s, "ArrowDown" )) return FTR_KEY_DOWN;
	if (0 == strcmp(s, "ArrowUp"   )) return FTR_KEY_UP;
	if (0 == strcmp(s, "ArrowLeft" )) return FTR_KEY_LEFT;
	if (0 == strcmp(s, "ArrowRight")) return FTR_KEY_RIGHT;
	if (0 == strcmp(s, "Delete"    )) return FTR_KEY_DEL;
	if (0 == strcmp(s, "Backspace" )) return FTR_KEY_BS;
	if (0 == strcmp(s, "Enter"     )) return '\n';
	if (0 == strcmp(s, "Shift"     )) return -1;
	if (0 == strcmp(s, "Control"   )) return -1;
	if (0 == strcmp(s, "Alt"       )) return -1;
	return 0;
}

//static int mod_from_emev(const EmscriptenKeyboardEvent *e)
//{
//	int m = 0;
//	if (e->shiftKey) m |= FTR_MASK_SHIFT;
//	if (e->ctrlKey ) m |= FTR_MASK_CONTROL;
//	if (e->altKey  ) m |= FTR_MASK_MOD1;
//	return m;
//}

static int mod_combine(int shift, int ctrl, int alt)
{
	int m = 0;
	if (shift) m |= FTR_MASK_SHIFT;
	if (ctrl ) m |= FTR_MASK_CONTROL;
	if (alt  ) m |= FTR_MASK_MOD1;
	return m;
}

// keyboard callback
EM_BOOL em_key(int type, const EmscriptenKeyboardEvent *e, void *usr)
{
	//printf("EM key event: type=%d key=%s code=%s 3us=%lu %lu %lu\n",
	//		type, e->key, e->code,
	//		e->charCode, e->keyCode, e->which
	//      );
	//fflush(stdout);

	struct _FTR *f = ftr_emscripten_global_state;

	f->last_m = mod_combine(e->shiftKey, e->ctrlKey, e->altKey);

	if (f->handle_key && type == EMSCRIPTEN_EVENT_KEYDOWN)
	{
		int k = key_from_emkey(e->key);
		if (k <= 0) return 1;
		if (e->shiftKey) k = toupper(k);
		f->handle_key((void*)f, k, f->last_m, f->last_x, f->last_y);
	}

	return 1;
}

static int button_from_emev(const EmscriptenMouseEvent *e)
{
	int b = 0;
	if (e->button == 0) b = FTR_BUTTON_LEFT;
	if (e->button == 1) b = FTR_BUTTON_MIDDLE;
	if (e->button == 2) b = FTR_BUTTON_RIGHT;
	return b;
}


// mouse callback
EM_BOOL em_mouse(int type, const EmscriptenMouseEvent *e, void *usr)
{
	printf("EM mouse event : type=%d x=%ld y=%ld button=%d buttons=%d",
			type, e->clientX, e->clientY, e->button, e->buttons);
	printf("\tmovementXY=%ld %ld targetXY=%ld %ld\n",
			e->movementX, e->movementY, e->targetX, e->targetY);
	fflush(stdout);

	struct _FTR *f = ftr_emscripten_global_state;

	f->last_x = e->targetX;
	f->last_y = e->targetY;

	if (f->handle_button && (type==EMSCRIPTEN_EVENT_MOUSEUP
				|| type==EMSCRIPTEN_EVENT_MOUSEDOWN)
			)
	{
		int b = button_from_emev(e);
		if (type == EMSCRIPTEN_EVENT_MOUSEUP)
			b *= -1;
		int m = mod_combine(e->shiftKey, e->ctrlKey, e->altKey);
		f->last_m = m;
		f->handle_button((void*)f, b, m, e->targetX, e->targetY);
	}

	if (f->handle_motion && type==EMSCRIPTEN_EVENT_MOUSEMOVE)
	{
		//int m = mod_combine(e->shiftKey, e->ctrlKey, e->altKey);
		//f->last_m = m;
		// for motion events, the modifier encodes buttons position
		int m = 0;
		if (e->buttons & 1) m |= FTR_BUTTON_LEFT;
		if (e->buttons & 2) m |= FTR_BUTTON_RIGHT;
		if (e->buttons & 4) m |= FTR_BUTTON_MIDDLE;
		f->handle_motion((void*)f, 0, m, e->targetX, e->targetY);
	}

	return 1;
}

// scroll callback (for some reason, separate from mouse)
EM_BOOL em_wheel(int type, const EmscriptenWheelEvent *e, void *usr)
{
	//printf("EM wheel event type=%d: deltaX=%f deltaY=%f deltaZ=%f deltaMode=%lu\n", type, e->deltaX, e->deltaY, e->deltaZ, e->deltaMode);
	//fflush(stdout);

	struct _FTR *f = ftr_emscripten_global_state;

	if (f->handle_button && type==EMSCRIPTEN_EVENT_WHEEL)
	{
		int b = e->deltaY > 0 ? FTR_BUTTON_DOWN : FTR_BUTTON_UP;
		f->handle_button((void*)f, b, f->last_m, f->last_x, f->last_y);
	}
	return 1;
}



// internal function to setup the global state
static void setup_emscripten_environment(void)
{
	struct _FTR *f = ftr_emscripten_global_state;

	emscripten_set_canvas_element_size("#canvas", f->w, f->h);
	emscripten_set_mousedown_callback("#canvas", NULL, 1, em_mouse);
	emscripten_set_mousemove_callback("#canvas", NULL, 1, em_mouse);
	emscripten_set_mouseup_callback("#canvas", NULL, 1, em_mouse);
	emscripten_set_wheel_callback("#canvas", NULL, 1, em_wheel);
	const char *target = EMSCRIPTEN_EVENT_TARGET_DOCUMENT;
	emscripten_set_keydown_callback(target, NULL, 1, em_key);
	emscripten_set_keyup_callback(target, NULL, 1, em_key);
}







void ftr_change_title(struct FTR *ff, char *s)
{
	(void)ff;
	(void)s;
}



// ftr_new_window_with_image_uint8_rgb {{{2
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	struct _FTR f[1];
	printf("asking ems window of size %d x %d\n", w, h);

	f->w = f->max_w = w;
	f->h = f->max_h = h;
	f->rgb = malloc(w * h * 3);
	f->rgba = malloc(w * h * 4);
	for (int i = 0; i < 3*w*h; i++)
		f->rgb[i] = x ? x[i] : 0;

	f->handle_key = NULL;//ftr_handler_exit_on_ESC;
	f->handle_button = NULL;
	f->handle_motion = NULL;
	f->handle_expose = NULL;
	f->handle_resize = NULL;
	f->handle_idle = NULL;
	f->handle_idle_toggled = NULL;
	f->stop_loop = 0;

	f->last_m = 0;
	f->last_x = -1;
	f->last_y = -1;


	return *(struct FTR *)f;
	//// avoid stupid type-punning warning:
	//struct FTR F[1];
	//memcpy(F, f, sizeof f);
	//return *F;
}


// ftr_close {{{2
void ftr_close(struct FTR *ff)
{
	fprintf(stderr, "FTR CLOSE!\n");
}

// ftr_loop_run {{{2
int ftr_loop_run(struct FTR *ff)
{
	fprintf(stdout, "going to start loop\n");
	struct _FTR *f = (void*)ff;

	ftr_emscripten_global_state = f;

	setup_emscripten_environment();

	fprintf(stdout, "entering loop run!\n");
	//if (f->handle_expose)
	//	f->handle_expose(ff, 0, 0, 0, 0);
	//ftr_dump_rgb_to_rgba(f);


	// the loop
	emscripten_set_main_loop(em_render, 0, 1);

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
