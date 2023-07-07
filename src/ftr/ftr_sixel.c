#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
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

// auxiliary stream interface, used mainly for sixels
struct bytestream {
	int n, ntop;
	uint8_t *t;
};

static void bytestream_init(struct bytestream *s)
{
	s->ntop = 1024;
	s->n = 0;
	s->t = malloc(s->ntop);
}

static void bs_putchar(struct bytestream *s, uint8_t x)
{
	if (s->n >= s->ntop)
	{
		s->ntop *= 2;
		s->t = realloc(s->t, s->ntop);
	}
	s->t[s->n++] = x;
}

static void bs_puts(struct bytestream *s, char *x)
{
	while (*x)
		bs_putchar(s, *x++);
}

static int bs_printf(struct bytestream *s, char *fmt, ...)
{
	va_list argp;
	char buf[0x1000];
	va_start(argp, fmt);
	int r = vsnprintf(buf, 0x1000, fmt, argp);
	bs_puts(s, buf);
	va_end(argp);
	return r;
}

static void bytestream_free(struct bytestream *s)
{
	free(s->t);
}

static int sidx(uint8_t *rgb) // sixel index identifier
{
	int r = 0;
	r += (rgb[0] >> 5) << 5;
	r += (rgb[1] >> 5) << 2;
	r += rgb[2] >> 6;
	return r;
	//return 32*(rgb[0]/32) + 4*(rgb[1]/32) + rgb[2]/64;
}

static void dump_sixels_to_bytestream_rgb3(
		struct bytestream *out,
		uint8_t *x, int w, int h)
{
	bs_puts(out, "\033Pq\n");
	for (int i = 0; i < 0x100; i++)
		bs_printf(out, "#%d;2;%d;%d;%d", i,
				(int)(14.2857*(i/32)),
				(int)(14.2857*((i/4)%8)),
				(int)(33.3333*(i%4)));
	for (int j = 0; j < h/6; j++)
	{
		int m[0x100] = {0}, c = 0;
		for (int i = 0; i < 6*w; i++)
		{
			int k = sidx(x+3*(6*j*w+i));
			if (!m[k]) c += 1;
			m[k] += 1;
		}
		for (int k = 0; k < 0x100; k++)
		if (m[k])
		{
			int b[w], r[w], R[w], n = 0;
			c -= 1;
			for (int i = 0; i < w; i++)
			{
				int s = 0;
				for (int l = 5; l >= 0; l--)
					s = 2*s + (k==sidx(x+3*((6*j+l)*w+i)));
				b[i] = s + 63;
			}
			for (int i = 0; i < w; i++)
				R[i] = 1;
			r[0] = *b;
			for (int i = 1; i < w; i++)
				if (b[i] == r[n])
					R[n] += 1;
				else
					r[++n] = b[i];
			bs_printf(out, "#%d", k);
			for (int i = 0; i <= n; i++)
				if (R[n] < 3)
					for (int l = 0; l < R[i]; l++)
						bs_putchar(out, r[i]);
				else
					bs_printf(out, "!%d%c", R[i], r[i]);
			bs_puts(out, c ? "$\n" : "-\n");
		}
	}
	bs_puts(out, "\033\\");
}

static void dump_sixels_to_bytestream_gray2(
		struct bytestream *out,
		uint8_t *x, int w, int h)
{
	int Q = (1<<2); // quantization over [0..255]
	bs_printf(out, "\033Pq\n");
	for (int i = 0; i < 0x100/Q; i++)
		bs_printf(out, "#%d;2;%d;%d;%d",
			i, (int)(Q*.39*i), (int)(Q*.39*i), (int)(Q*.39*i));
	for (int j = 0; j < h/6; j++) {
		int m[0x100] = {0}, c = 0;
		for (int i = 0; i < 6*w; i++) {
			int k = x[6*j*w+i]/Q;
			if (!m[k]) c += 1;
			m[k] += 1;
		}
		for (int k = 0; k < 0x100/Q; k++)
		if (m[k]) {
			int b[w], r[w], R[w], idx = 0;
			c -= 1;
			for (int i = 0; i < w; i++) {
				b[i] = 0;
				for (int l = 5; l >= 0; l--)
					b[i] = 2*b[i] + (k == x[(6*j+l)*w+i]/Q);
				b[i] += 63;
			}
			for (int i = 0; i < w; i++) R[i] = 1;
			r[0] = *b;
			for (int i = 1; i < w; i++)
				if (b[i] == r[idx]) R[idx] += 1;
				else r[++idx] = b[i];
			bs_printf(out, "#%d", k);
			for (int i = 0; i <= idx; i++)
				if (R[idx] < 3)
					for (int l = 0; l < R[i]; l++)
						bs_printf(out, "%c", r[i]);
				else
					bs_printf(out, "!%d%c", R[i], r[i]);
			bs_printf(out, c ? "$\n" : "-\n");
		}
	}
	bs_printf(out, "\033\\");
}

//static void dump_sixels_to_stdout_rgb3(uint8_t *x, int w, int h)
//{
//	struct bytestream s[1];
//	bytestream_init(s);
//	dump_sixels_to_bytestream_rgb3(s, x, w, h);
//	for (int i = 0; i < s->n; i++)
//		putchar(s->t[i]);
//	bytestream_free(s);
//	//printf("\ePq\n");
//	//for (int i = 0; i < 0x100; i++)
//	//	printf("#%d;2;%d;%d;%d", i,(int)(14.2857*(i/32)),
//	//			(int)(14.2857*((i/4)%8)), (int)(33.3333*(i%4)));
//	//for (int j = 0; j < h/6; j++)
//	//{
//	//	int m[0x100] = {0}, c = 0;
//	//	for (int i = 0; i < 6*w; i++)
//	//	{
//	//		int k = sidx(x+3*(6*j*w+i));
//	//		if (!m[k]) c += 1;
//	//		m[k] += 1;
//	//	}
//	//	for (int k = 0; k < 0x100; k++)
//	//	if (m[k])
//	//	{
//	//		int b[w], r[w], R[w], n = 0;
//	//		c -= 1;
//	//		for (int i = 0; i < w; i++)
//	//		{
//	//			int s = 0;
//	//			for (int l = 5; l >= 0; l--)
//	//				s = 2*s + (k==sidx(x+3*((6*j+l)*w+i)));
//	//			b[i] = s + 63;
//	//		}
//	//		for (int i = 0; i < w; i++)
//	//			R[i] = 1;
//	//		r[0] = *b;
//	//		for (int i = 1; i < w; i++)
//	//			if (b[i] == r[n])
//	//				R[n] += 1;
//	//			else
//	//				r[++n] = b[i];
//	//		printf("#%d", k);
//	//		for (int i = 0; i <= n; i++)
//	//			if (R[n] < 3)
//	//				for (int k = 0; k < R[i]; k++)
//	//					printf("%c", r[i]);
//	//			else
//	//				printf("!%d%c", R[i], r[i]);
//	//		printf(c ? "$\n" : "-\n");
//	//	}
//	//}
//	//printf("\e\\");
//}

//static void dump_sixels_to_stdout_gray2(uint8_t *x, int w, int h)
//{
//	int Q = (1<<2); // quantization over [0..255]
//	printf("\ePq\n");
//	for (int i = 0; i < 0x100/Q; i++)
//		printf("#%d;2;%d;%d;%d",
//			i, (int)(Q*.39*i), (int)(Q*.39*i), (int)(Q*.39*i));
//	for (int j = 0; j < h/6; j++) {
//		int m[0x100] = {0}, c = 0;
//		for (int i = 0; i < 6*w; i++) {
//			int k = x[6*j*w+i]/Q;
//			if (!m[k]) c += 1;
//			m[k] += 1;
//		}
//		for (int k = 0; k < 0x100/Q; k++)
//		if (m[k]) {
//			int b[w], r[w], R[w], idx = 0;
//			c -= 1;
//			for (int i = 0; i < w; i++) {
//				b[i] = 0;
//				for (int l = 5; l >= 0; l--)
//					b[i] = 2*b[i] + (k == x[(6*j+l)*w+i]/Q);
//				b[i] += 63;
//			}
//			for (int i = 0; i < w; i++) R[i] = 1;
//			r[0] = *b;
//			for (int i = 1; i < w; i++)
//				if (b[i] == r[idx]) R[idx] += 1;
//				else r[++idx] = b[i];
//			printf("#%d", k);
//			for (int i = 0; i <= idx; i++)
//				if (R[idx] < 3)
//					for (int k = 0; k < R[i]; k++)
//						printf("%c", r[i]);
//				else
//					printf("!%d%c", R[i], r[i]);
//			printf(c ? "$\n" : "-\n");
//		}
//	}
//	printf("\e\\");
//}

static void penetrate_screen(struct bytestream *z, struct bytestream *s)
{
	//bs_puts(z, "PENETRATE SCREEN\n");
	bs_puts(z, "\x1bP");
	for (int i = 0; i < s->n; i++)
	{
		int c = s->t[i];
		     if (c == '\x90') bs_puts(z, "\x1bP");
		else if (c == '\x9c') bs_puts(z, "\x1b\x1b\\\x1bP\\");
		else if (c == '\x1b' && i+1 < s->n && s->t[i] == '\\')
		{
			bs_puts(z, "\x1b\x1b\\\x1bP\\");
			i += 1;
		}
		else bs_putchar(z, c);
	}
	bs_puts(z, "\x1b\\");
}

static void dump_sixels_to_stdout_uint8(uint8_t *x, int w, int h, int pd)
{

	struct bytestream s[1];
	bytestream_init(s);
	if (pd == 3) dump_sixels_to_bytestream_rgb3(s, x, w, h);
	if (pd == 1) dump_sixels_to_bytestream_gray2(s, x, w, h);
	struct bytestream *S = s;
	bool screen = false;//strstr(xgetenv("TERM"), "screen");
	struct bytestream z[1];
	if (screen) bytestream_init(z);
	if (screen) penetrate_screen(S = z, s);
	for (int i = 0; i < S->n; i++)
		putchar(S->t[i]);
	bytestream_free(s);
	if (screen) bytestream_free(z);
}

//static void dump_sixels_to_stdout(struct iio_image *x)
//{
//	//if (x->type != IIO_TYPE_UINT8)
//	{
//		void *old_data = x->data;
//		int ss = iio_image_sample_size(x);
//		int nsamp = iio_image_number_of_samples(x);
//		x->data = malloc(nsamp*ss);
//		memcpy(x->data, old_data, nsamp*ss);
//		iio_convert_samples(x, IIO_TYPE_UINT8);
//		dump_sixels_to_stdout_uint8(x->data, x->sizes[0], x->sizes[1],
//				x->pixel_dimension);
//		//if (x->pixel_dimension==3)
//		//dump_sixels_to_stdout_rgb3(x->data, x->sizes[0], x->sizes[1]);
//		//else if (x->pixel_dimension==1)
//		//dump_sixels_to_stdout_gray2(x->data,x->sizes[0],x->sizes[1]);
//		free(x->data);
//		x->data = old_data;
//	}
//}





static void idump_six(uint8_t *x, int w, int h)
{
	// dump an image at the top of the screen using ANSI sequences
	// NOTE: 0x1b == 033 == 27 == ESC

	printf("\x1b[1;1H");
	dump_sixels_to_stdout_uint8(x, w, h, 3);
}

static void ftr_term_dump(struct _FTR *f)
{
	printf("dump %d %d:\n", f->w, f->h);
	idump_six(f->rgb, f->w, f->h);
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
	w = 800;
	h = 600;

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

	//return *(struct FTR *)f;
	// avoid stupid type-punning warning:
	struct FTR F[1];
	memcpy(F, f, sizeof f);
	return *F;
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

void ftr_x11_force_redraw(struct FTR *ff)
{
	// just ignore it
	(void)ff;
}

#include "ftr_common_inc.c"

// vim:set foldmethod=marker:
