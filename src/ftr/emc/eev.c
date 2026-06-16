#include <emscripten.h>
#include <emscripten/html5.h>
#include <stdio.h>

// Simple framebuffer
static const int global_W = 800;
static const int global_H = 600;
static unsigned char global_pixels[global_W * global_H * 4];
static int global_dirty = 1;

void green_background()
{
	for (int i = 0; i < global_W * global_H; i++) {
		global_pixels[i * 4 + 0] = 0;   // R
		global_pixels[i * 4 + 1] = 155; // G
		global_pixels[i * 4 + 2] = 0;   // B
		global_pixels[i * 4 + 3] = 255; // A
	}
}

// no arguments ! (we need global variables)
void render(void)
{
	if (!global_dirty) return;
	global_dirty = 0;

	EM_ASM({
		const c = Module.canvas.getContext('2d');
		const p = $0; // pixels
		const w = $1; // W
		const h = $2; // H
		const v = new Uint8ClampedArray(Module.HEAPU8.buffer, p, w*h*4);
		const i = new ImageData(v, w, h);
		c.putImageData(i, 0, 0);
		},
	global_pixels, global_W, global_H);
}

// mouse callback
EM_BOOL mouse_cb(int eventType, const EmscriptenMouseEvent *e, void *usr)
{
	printf("mouse event : type=%d x=%ld y=%ld button=%d\n",
			eventType, e->clientX, e->clientY, e->button);
	fflush(stdout);
	return 1;
}

// scroll callback (for some reason, separate from mouse)
EM_BOOL wheel_cb(int eventType, const EmscriptenWheelEvent *e, void *usr)
{
	printf("wheel event : deltaX=%f deltaY=%f deltaZ=%f deltaMode=%lu\n",
			e->deltaX, e->deltaY, e->deltaZ, e->deltaMode);

	fflush(stdout);
	return 1;
}

// keyboard callback
EM_BOOL key_cb(int eventType, const EmscriptenKeyboardEvent *e, void *usr)
{
	printf("Key event: type=%d key=%s code=%s\n",
			eventType, e->key, e->code);
	fflush(stdout);
	return 1;
}

int main() {

	// set canvas size
	emscripten_set_canvas_element_size("#canvas", global_W, global_H);

	// paint into local canvas
	green_background();

	// mouse
	emscripten_set_mousedown_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_mousemove_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_mouseup_callback("#canvas", NULL, 1, mouse_cb);
	emscripten_set_wheel_callback( "#canvas", NULL, 1, wheel_cb);

	// keyboard
	const char *target = EMSCRIPTEN_EVENT_TARGET_DOCUMENT; // #canvas has no focus
	emscripten_set_keydown_callback(target, NULL, 1, key_cb);
	emscripten_set_keyup_callback(target, NULL, 1, key_cb);

	// render loop (moves the local canvas to the browser)
	emscripten_set_main_loop(render, 0, 1);

	return 0;
}
