// str: simpler terminal
//
//     That the concept of 'pty' still exists in the year 2009 is quite fucking
//     amazing. I'm surprised we don't carry punchcards around anymore.
//     --uriel
//
//
// This is not a terminal emulator.
// It does not have a concept of "tty" or "session", or other bullshit.
// This is just a program that shows the output of another program as
// characters drawn on a rectangular window.
// The goal of this program is to be as simple as possible and fulfill the
// required tast; there is no concern for efficiency.
//

#include "ftr.c"
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "xfont9x15.c"

// a terminal is just a 80x25 matrix of characters, with some options
struct terminal {
	// essential data
	int w, h; // size in characters (always 80x25)
	int *letters;    // unicode points

	// ancillary
	int cursorx, cursory;
	int kerning;
	int spacing;
	struct bitmap_font font[1];
	int *attributes; // colors and whatnot
};

// affect the terminal state from the given string
void term_puts(struct terminal *t, char *s)
{
	t->letters[0] = 'a';
	t->letters[1] = 'b';
	t->letters[2] = 'c';
	while (1)
	{
		int c = *s++;
		if (!c) break;
		if (c == '\n') {
		///etc
		}
	}
}

void term_add_char_under_cursor(struct terminal *t, int c)
{
	if (t->cursorx >= 0 && t->cursorx < t->w)
	if (t->cursory >= 0 && t->cursory < t->h)
		t->letters[t->cursory * t->w + t->cursorx] = c;
}

void term_new_line(struct terminal *t)
{
	t->cursory += 1;
	t->cursorx = 0;
}


void term_advance_cursor(struct terminal *t)
{
	if (t->cursorx < t->w -1)
		t->cursorx += 1;
	else
		term_new_line(t);
}

static void
put_char_in_rgb_image(uint8_t *rgb, int w, int h, struct bitmap_font *f,
		int posx, int posy, uint8_t *fg, uint8_t *bg, int c)
{
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

	// clear rgb buffer
	for (int i = 0; i < 3 * w * h; i++)
		rgb[i] = 0;

	// draw character bitmaps
	for (int j = 0; j < t->h; j++)
	for (int i = 0; i < t->w; i++)
	{
		int c = t->letters[j * t->w + i];
		int a = t->attributes[j * t->w + i];
		uint8_t fg[3] = {255, 255, 255};
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

#include <ctype.h>
static void term_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	fprintf(stderr, "TERM_KEY_HANDLER  %d '%c' (%d) at %d %d\n",
			k, isalpha(k)?k:' ', m, x, y);
	if  (k == '\033')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	// todo: do not echo the keys directly, sent to the terminal input
	// stream, just like the output of the program

	struct terminal *t = f->userdata;
	if (k == FTR_KEY_LEFT ) t->cursorx -= 1;
	if (k == FTR_KEY_RIGHT) t->cursorx += 1;
	if (k == FTR_KEY_UP   ) t->cursory -= 1;
	if (k == FTR_KEY_DOWN ) t->cursory += 1;

	if (isprint(k) || k==' ') {
		if (m & FTR_MASK_SHIFT && isalpha(k))
			k = toupper(k);
		term_add_char_under_cursor(t, k);
		term_advance_cursor(t);
	}

	if (k == '\b' && t->cursorx > 0)
	{
		t->cursorx -= 1;
		term_add_char_under_cursor(t, ' ');
	}

	if (k == '\n')
		term_new_line(t);

	f->changed = 1;
}

static void term_exposer(struct FTR *f, int b, int m, int x, int y)
{
	fprintf(stderr, "term expose\n");
	term_bitmap(f->rgb, f->w, f->h, (struct terminal *)f->userdata);
	f->changed = 1;
}

int main()
{
	struct terminal t[1];
	t->w = 80;
	t->h = 25;
	t->kerning = 0;
	t->spacing = 0;
	t->font[0] = *xfont9x15;
	t->font[0] = reformat_font(t->font[0], UNPACKED);
	t->letters = malloc(sizeof(int) * t->w * t->h);
	t->attributes = malloc(sizeof(int) * t->w * t->h);
	t->cursorx = t->cursory = 0;
	int w = (t->w + t->kerning) * t->font->width;
	int h = (t->h + t->spacing) * t->font->height;

	term_puts(t, "");

	struct FTR f = ftr_new_window(w, h);
	f.userdata = t;
	f.changed = 1;

	ftr_set_handler(&f, "expose", term_exposer);
	ftr_set_handler(&f, "key", term_key_handler);
	return ftr_loop_run(&f);
}
