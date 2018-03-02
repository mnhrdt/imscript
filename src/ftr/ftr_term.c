#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>

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

};

// Check that _FTR can fit inside a FTR
// (if this line fails, increase the padding at the end of struct FTR on ftr.h)
typedef char check_FTR_size[sizeof(struct _FTR)<=sizeof(struct FTR)?1:-1];


static void idump(uint8_t *x, int w, int h)
{
	// dump an image at the top of the screen using ANSI sequences
	// NOTE: 0x1b == 033 == 27 == ESC

	printf("\x1b[1;1H");
	for (int j = 0; j < h; j += 2)
	{
		for (int i = 0; i < w; i++)
		{
			uint8_t *a = x + 3 * (w*(j+0) + i);
			uint8_t *b = x + 3 * (w*(j+1) + i);
			printf("\x1b[38;2;%d;%d;%dm\x1b[48;2;%d;%d;%dm▀",
					a[0], a[1], a[2], b[0], b[1], b[2]);
		}
		printf("\x1b[0m\n");
	}
}

static void ftr_term_dump(struct _FTR *f)
{
	printf("dump %d %d:\n", f->w, f->h);
	idump(f->rgb, f->w, f->h);
}

static void disable_canonical_and_echo_modes(void)
{
	struct termios t[1];
	tcgetattr(0, t);
	t->c_lflag &= ~ICANON;
	t->c_lflag &= ~ECHO;
	tcsetattr(0, TCSANOW, t);
}

static void enable_canonical_and_echo_modes(void)
{
	struct termios t[1];
	tcgetattr(0, t);
	t->c_lflag |= ICANON;
	t->c_lflag |= ECHO;
	tcsetattr(0, TCSANOW, t);
}

void ftr_change_title(struct FTR *ff, char *s)
{
}

#include "smapa.h"
SMART_PARAMETER(COLUMNS,80)
SMART_PARAMETER(LINES,25)

// ftr_new_window_with_image_uint8_rgb {{{2
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	//if (w != 80 || h != 25)
	//	exit(fprintf(stderr, "this is not a proper terminal!\n"));
	w = COLUMNS();
	h = 2*(LINES() - 7) - 4;

	struct _FTR f[1];

	f->w = w;
	f->h = h;
	f->max_w = 2000;
	f->max_h = 2000;
	f->rgb = malloc(f->max_w * f->max_h * 3);
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

	disable_canonical_and_echo_modes();
	ftr_term_dump(f);

	return *(struct FTR *)f;
}

// ftr_close {{{2
void ftr_close(struct FTR *ff)
{
	fprintf(stderr, "FTR CLOSE!\n");
	enable_canonical_and_echo_modes();

	struct _FTR *f = (void*)ff;
	if (f->rgb) free(f->rgb);
}

// ftr_loop_run {{{2
int ftr_loop_run(struct FTR *ff)
{
	//fprintf(stderr, "going to start loop\n");
	struct _FTR *f = (void*)ff;

	fprintf(stderr, "entering loop run!\n");
	if (f->handle_expose)
		f->handle_expose(ff, 0, 0, 0, 0);
	ftr_term_dump(f);

	// only events are obtained by GETCHAR
	// they always redraw everything

	int buf_n = 0x100, buf[buf_n], buf_i = 0;
	bool inside_code = false;
	while (!f->stop_loop)
	{
		int c = getchar();
		if (c == EOF) break;
		fprintf(stderr, "RAW key %d '%c' (in=%s,idx=%d)\n",
				c, isprint(c)?c:' ',
				inside_code?"yes":"no", buf_i);

		if (inside_code) {
			if (buf_i == 0 && c == '[') { // first char of code
				buf[buf_i++] = c;
				continue;
			}
			else if (buf_i == 1 && buf[0] == '[' && isupper(c))
			{
				buf[buf_i++] = c;
			}
		} else {
			if (c == 27) // ESC
			{
				inside_code = true;
				continue;
			}
		}

		if (inside_code)
		{
			assert(buf_i == 2);
			if (buf[0] == '[' && buf[1] == 'A') c = FTR_KEY_UP;
			if (buf[0] == '[' && buf[1] == 'B') c = FTR_KEY_DOWN;
			if (buf[0] == '[' && buf[1] == 'C') c = FTR_KEY_RIGHT;
			if (buf[0] == '[' && buf[1] == 'D') c = FTR_KEY_LEFT;
			buf_i = 0;
			inside_code = false;
		}

		if (c == 127) c = '\b';

		fprintf(stderr, "trans key %d '%c'\n", c, isprint(c)?c:' ');
		if (f->handle_key)
			f->handle_key(ff, c, 0, f->w/2, f->h/2);
		if (f->handle_expose)
			f->handle_expose(ff, 0, 0, 0, 0);
		ftr_term_dump(f);
	}

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

#include "ftr_common_inc.c"

// vim:set foldmethod=marker:
