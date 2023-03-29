// str: simpler terminal
//
//     That the concept of 'pty' still exists in the year 2009 is quite fucking
//     amazing. I'm surprised we don't carry punchcards around anymore.
//     --uriel
//
//
// This is not a terminal emulator.
//
// It does not have a concept of "tty" or "session", or other bullshit.
// This is just a program that shows the output of another program as
// characters drawn on a rectangular window.
// It does not use X fonts, only a simple function to dump a bitmap into a
// window.
// The goal of this program is to be as simple as possible and fulfill the
// required task; there is no concern for efficiency.
// The only (inevitable) legacy is that old vt100 control sequences are used.
// This is to avoid having to define an essentially arbitrary terminfo.
//
// This should come as fresh air in the world of bloated programs such as
// xterm, rxvt, or st.  Even the leanest of these programs (st) relies on huge
// dependences like libfreetype and obnoxious concepts such as pty.
//

#include <stdio.h> // fprintf, popen
#include <ctype.h> // isalpha, isprint

#include "ftr.h"
#define OMIT_MAIN_FONTU
#include "fontu.c" // todo: cherry-pick the required fontu functions
#include "fonts/xfont_10x20.c"

// a terminal is just a 80x25 matrix of characters, with some options
struct terminal {
	// essential data
	int w, h;        // size in characters (always 80x25)
	int *letters;    // table of W*H unicode points

	// ancillary
	int cursorx, cursory;
	int kerning;
	int spacing;
	struct bitmap_font font[1];
	int *attributes; // colors and whatnot (table of W*H)

	// output of the program running inside the terminal at this moment
	FILE *stream;
	FILE *streams[2];

	// program management
	int program;
};

// Overall structure:
// 1. a program P is run
// 2. the stdout and stderr of P are dumped into the state machine
// 3. keyboard events are dumped into the stdin of P
// 4. unless "canonicalized", keyboard events are also dumped into the state
// machine


// actions

void term_action_put_char_under_cursor(struct terminal *t, int c)
{
	if (t->cursorx >= 0 && t->cursorx < t->w)
	if (t->cursory >= 0 && t->cursory < t->h)
		t->letters[t->cursory * t->w + t->cursorx] = c;
}

void term_action_new_line(struct terminal *t)
{
	t->cursory += 1;
	t->cursorx = 0;
}

void term_action_cursor_dxy(struct terminal *t, int dx, int dy)
{
	t->cursorx += dx;
	t->cursory += dy;
	fprintf(stderr, "CURSOR dxy (%d %d), now is (%d %d)\n", dx, dy,
			t->cursorx, t->cursory);
}

void term_action_advance_cursor(struct terminal *t)
{
	if (t->cursorx < t->w -1)
		t->cursorx += 1;
	else
		term_action_new_line(t);
}

// affect the terminal state from the given string
//
// (The terminal is a finite state machine.  This is the only function
// that can be used to edit the state.  All changes go thru calls to
// functions of the form "term_action_*".)
void term_puts(struct terminal *t, char *s)
{
	fprintf(stderr, "term puts s=\"%s\"\n", s);
	while (1)
	{
		int c = *s++;
		if (!c) break;
		else if (c == '\n')
			term_action_new_line(t);
		else if (c == 0x1b) {
			c = *s++;
			if (!c) break;
			if (c == 'A') term_action_cursor_dxy(t,  0, -1);
			if (c == 'B') term_action_cursor_dxy(t,  0,  1);
			if (c == 'C') term_action_cursor_dxy(t,  1,  0);
			if (c == 'D') term_action_cursor_dxy(t, -1,  0);
		}
		else {
			term_action_put_char_under_cursor(t, c);
			term_action_advance_cursor(t);
		}
	}
}

static void
put_char_in_rgb_image(uint8_t *rgb, int w, int h, struct bitmap_font *f,
		int posx, int posy, uint8_t *fg, uint8_t *bg, int c)
{
	if (!c) return;
	if (c > 0 && c < f->number_of_glyphs)
		for (int i = 0; i < f->width; i++)
		for (int j = 0; j < f->height; j++)
		{
			int ii = posx + i;
			int jj = posy + j;
			if (get_font_bit(f, c, i, j))
				put_pixel_rgb(rgb, w, h, ii, jj, fg);
			else
				put_pixel_rgb(rgb, w, h, ii, jj, bg);
		}
}

// dump the terminal state into a rgb bitmap
void term_bitmap(uint8_t *rgb, int w, int h, struct terminal *t)
{
	// compute grid offsets
	int ox = t->font->width  + t->kerning;
	int oy = t->font->height + t->spacing;

	// size consistency check
	if (w != ox * t->w || h != oy * t->h)
		fprintf(stderr, "WARNING: (%d,%d) != (%d,%d)\n",
				w, h, ox*t->w, oy*t->h);

	// clear rgb buffer
	for (int i = 0; i < 3 * w * h; i++)
		rgb[i] = 0;

	// draw character bitmaps
	for (int j = 0; j < t->h; j++)
	for (int i = 0; i < t->w; i++)
	{
		int c = t->letters[j * t->w + i];
		int a = t->attributes[j * t->w + i];
		uint8_t fg[3] = {205, 205, 205};
		uint8_t bg[3] = {0, 0, 0};
		put_char_in_rgb_image(rgb,w,h, t->font, i*ox,j*oy, fg,bg, c);
	}

	// draw cursor
	if (t->cursorx >= 0 && t->cursorx < t->w)
	if (t->cursory >= 0 && t->cursory < t->h)
	{
		int c0x = t->cursorx * ox;
		int c0y = t->cursory * oy;
		uint8_t cc[3] = {200, 100, 0};
		for (int i = 0; i < ox; i++) {
			put_pixel_rgb(rgb,w,h, c0x + i, c0y     , cc);
			put_pixel_rgb(rgb,w,h, c0x + i, c0y + oy, cc);
		}
		for (int j = 0; j < oy; j++) {
			put_pixel_rgb(rgb,w,h, c0x     , c0y + j, cc);
			put_pixel_rgb(rgb,w,h, c0x + ox, c0y + j, cc);
		}
	}
}


//static int bidirectional_pipe(FILE *outin[2], const char *program)
//{
//}
//
//static FILE *unidirectional_read_pipe(const char *program)
//{
//}

static FILE *my_popen_read(const char *s)
{
	int p[2]; // file descriptors of the pipe: p[0]=read p[1]=write
	if (pipe(p) != 0)
		return 0;

	int pid = fork();
	if (pid == 0) {          // CHILD
		close(p[0]);     // close read-end of the pipe
		dup2(p[1], 1);   // attach 1=stdout to write-end of the pipe
		close(p[1]);     // close the write-end of the pipe
		execl(s, "s", "caca", NULL);
		return 0;
	} else {                 // PARENT, pid = pid of child
		close(p[1]);     // close write-end of the pipe
		return fdopen(p[0], "r");
	}
}

static int my_popen_readwrite(FILE *rw[2], const char *s)
{
	int pr[2], pw[2];
	if (pipe(pr)) return 1;
	if (pipe(pw)) return 2;

	fprintf(stderr, "pr = %d %d\n", pr[0], pr[1]);
	fprintf(stderr, "pw = %d %d\n", pw[0], pw[1]);

	int pid = fork();
	fprintf(stderr, "pid = %d\n", pid);
	if (pid == 0) { // CHILD
		close(pr[1]);
		close(pw[0]);
		dup2(pr[0], 0); // stdin
		dup2(pw[1], 1); // stdout
		fprintf(stderr, "i'm the child!\n");
		execl(s, "s", "caca", NULL);
		fprintf(stderr, "the child is done\n");
		_exit(127);
	}
	// PARENT, pid = pid of child
	close(pr[0]);
	close(pw[1]);
	fprintf(stderr, "i'm the parent!\n");
	fprintf(stderr, "rw[0] %d, rw[1] %d\n", pr[1], pw[0]);
	rw[0] = fdopen(pr[1], "w");
	rw[1] = fdopen(pw[0], "r");
	return 0;
}

static void term_action_run_dummy(struct terminal *t)
{
	if (t->program) return;
	t->program = 1;
	fprintf(stderr, "TERM ACTION RUN DUMMY\n");

	//t->stream  = popen("/bin/yes", "r");
	//t->stream  = my_popen_read("/bin/yes");
	//t->stream  = popen("./dummy2 write", "r");
	int r = my_popen_readwrite(t->streams, "./dummy");
	fprintf(stderr, "my_popen_readwrite = %d\n", r);
}

static void term_action_getchar(struct terminal *t)
{
	if (!t->program) return;
	fprintf(stderr, "term action getchar\n");

	fputc('p', t->streams[0]);
	fprintf(stderr, "   just gave it 'p'!\n");
	fflush(t->streams[0]); fflush(t->streams[1]);
	int c = fgetc(t->streams[1]);
	fflush(t->streams[0]); fflush(t->streams[1]);
	fprintf(stderr, "stream gave me c=%d '%c'\n", c, isprint(c)?c:' ');

	char s[2];
	s[0] = c;
	s[1] = 0;
	term_puts(t, s);
}

static char *global_table_of_keycodes_to_consolecodes[0x2000] = {
	[FTR_KEY_UP]    = "A",
	[FTR_KEY_DOWN]  = "B",
	[FTR_KEY_RIGHT] = "C",
	[FTR_KEY_LEFT]  = "D",
};

static void term_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "TERM_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
	if (k == '\033')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);
	if (k == 'x') {
		term_action_run_dummy(f->userdata);
		return;
	}
	if (k == ' ') {
		term_action_getchar(f->userdata);
		f->changed = 1;
		return;
	}


	if (k >= 0x2000) {
		fprintf(stderr, "WARNING: rejected key %d\n", k);
		return;
	}

	// 1. prepare a buffer representing this key
	char buf[0x10] = {0};

	// TODO: build a table of keys and strings (about 2000 entries)
	char *tk = global_table_of_keycodes_to_consolecodes[k];
	if (k < 0x2000 && tk && tk[0] == 0x1b)
		strcpy(buf, tk);

	if (isprint(k) || k==' ' || k=='\n')
	{
		if (m & FTR_MASK_SHIFT && isalpha(k))
			k = toupper(k);
		buf[0] = k;
	}


	// 2. dump this buffer into the state machine
	struct terminal *t = f->userdata;
	term_puts(t, buf);

	f->changed = 1;
}

static void term_exposer(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "term expose\n");
	// todo 1: expose only the changed part
	// todo 2: if scroll, copy without recomputing
	term_bitmap(f->rgb, f->w, f->h, (struct terminal *)f->userdata);
	f->changed = 1;
}

int main_srt()
{
	struct terminal t[1];
	t->w = 80;
	t->h = 25;
	t->kerning = 0;
	t->spacing = 0;
	t->font[0] = reformat_font(*xfont_10x20, UNPACKED);
	fprintf(stderr, "NUMBER_OF_GLYPHS = %d\n", t->font->number_of_glyphs);
	t->letters    = calloc(t->w * t->h, sizeof*t->letters);
	t->attributes = calloc(t->w * t->h, sizeof*t->attributes);
	t->cursorx = t->cursory = 0;
	int w = (t->w + t->kerning) * t->font->width;
	int h = (t->h + t->spacing) * t->font->height;

	t->program = 0;

	term_puts(t, "abc");

	struct FTR f = ftr_new_window(w, h);
	f.userdata = t;
	f.changed = 1;

	ftr_set_handler(&f, "expose", term_exposer);
	ftr_set_handler(&f, "key", term_key_handler);
	return ftr_loop_run(&f);
}

int main() { return main_srt(); }

// philosophy:
//
// 0. PIXEL IS KING: the terminal window is a grid of pixels
// 1. NO BLURRY FONTS: only monospace bitmap fonts are supported
// 2. NO SCALING: a pixel is a pixel.  If text is too small, use a larger font
// 3. NO GPU: the terminal runs in a single cpu thread
// 4. NO LIBRARIES: compile it with 0 dependencies and dump output to ppm files
// 5. NATIVE IMAGES: direct access to the terminal framebuffer via control codes
// 6. NATIVE VECTOR PRIMITIVES: for simple in-band plotting (segment, polygon)
// 7. MOCK XTERM: built silently on top of TERM=xterm-256color, no terminfo
// 8. LIMITED ENVELOPE: can run bash, vim, tmux, bpython and that's it
