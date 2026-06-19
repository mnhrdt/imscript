#include <ctype.h>
#include <stdio.h>

#include "ftr.h"
#include "random.c"

// bitmap fonts
#define OMIT_MAIN_FONTU
#include "fontu.c"
#include "fonts/xfonts_all.c"


struct bitmap_font font[1];

static void print_event_key(struct FTR *f, int k, int m, int x, int y)
{
	uint8_t red[3] = {255, 0, 0};
	uint8_t green[3] = {0, 255, 0};
	uint8_t dgreen[3] = {0, 155, 0};
	uint8_t dgray[3] = {100,100,100};
	uint8_t white[3] = {255, 255, 255};
	uint8_t dblue[3] = {0, 0, 227};
	uint8_t cyan[3] = {0, 255, 255};
	uint8_t black[3] = {0, 0, 0};
	uint8_t pink[3] = {255, 0, 255};
	char buf[0x100];
	snprintf(buf, 0x100, "event KEY     k=%d '%c'\nm=%d (%d %d)\n", k,
			(isprint(k))?k:' ', m, x, y);
	put_string_in_rgb_image(f->rgb, f->w, f->h,
		10, 10, green, black, 0, font, buf);

	//int r = 255*random_uniform();
	//int g = 255*random_uniform();
	//int b = 255*random_uniform();
	//printf("rgb = %d %d %d\n", r, g, b);
	//for (int i = 0; i < f->w * f->h; i++)
	//{
	//	f->rgb[3*i+0] = r;
	//	f->rgb[3*i+1] = g;
	//	f->rgb[3*i+2] = b;
	//}
	f->changed = 1;
}
static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && i < w && j >= 0 && j < h;
}
static void print_event_button(struct FTR *f, int k, int m, int x, int y)
{
	printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
}
static void print_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	//printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
	if (insideP(f->w, f->h, x, y))
	{
		f->rgb[3*(x+y*f->w)+0] = 255;
		f->rgb[3*(x+y*f->w)+1] = 255;
		f->rgb[3*(x+y*f->w)+2] = 255;
	}
	f->changed = 1;
}
static void expose(struct FTR *f, int ev_b, int ev_m, int ev_x, int ev_y)
{
	;
}

//#include "pickopt.c"
//int main_jmgs(int c, char *v[])
int main(void)//int c, char *v[])
{
	int w = 600;//atoi(pick_option(&c, &v, "w", "800"));
	int h = 600;//atoi(pick_option(&c, &v, "h", "800"));

	font[0] = reformat_font(*xfont_9x18B, UNPACKED);

	struct FTR f = ftr_new_window(600,600);
	for (int i = 0; i < 3 * f.w * f.h; i++)
		f.rgb[i] = 127;//0xa0*!(i%3);
	f.changed = 1;
	fprintf(stderr, "i'm here!\n");
	ftr_set_handler(&f, "key", print_event_key);
	ftr_set_handler(&f, "button", print_event_button);
	ftr_set_handler(&f, "motion", print_event_motion);
	ftr_set_handler(&f, "expose", expose);
	fprintf(stderr, "i'm still here!\n");
	ftr_loop_run(&f);
	return 0;
}
