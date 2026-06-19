#include <ctype.h>
#include <stdio.h>
#include "ftr.h"
#include "random.c"
static void print_event_key(struct FTR *f, int k, int m, int x, int y)
{
	printf("event KEY     k=%d '%c'\tm=%d (%d %d)\n", k,
			(isprint(k))?k:' ', m, x, y);
	int r = 255*random_uniform();
	int g = 255*random_uniform();
	int b = 255*random_uniform();
	printf("rgb = %d %d %d\n", r, g, b);
	for (int i = 0; i < f->w * f->h; i++)
	{
		f->rgb[3*i+0] = r;
		f->rgb[3*i+1] = g;
		f->rgb[3*i+2] = b;
	}
	f->changed = 1;
}
static void print_event_button(struct FTR *f, int k, int m, int x, int y)
{
	printf("event BUTTON  b=%d\tm=%d (%d %d)\n", k, m, x, y);
}
static void print_event_motion(struct FTR *f, int b, int m, int x, int y)
{
	printf("event MOTION  b=%d\tm=%d (%d %d)\n", b, m, x, y);
}
static void expose(struct FTR *f, int ev_b, int ev_m, int ev_x, int ev_y)
{
	;
}
int main(void)
{
	struct FTR f = ftr_new_window(600,600);
	for (int i = 0; i < 3 * f.w * f.h; i++)
		f.rgb[i] = 0xa0*!(i%3);
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
