#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include "ftr.h"
#include "cam.h"

struct kam_state {
	int dummy;
	struct camera_t *c;
};

static void kam_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	(void)unused_x; (void)unused_y;
	//fprintf(stderr, "\n\nexpose %d %d\n", b, m);
	struct kam_state *e = f->userdata;
	//fprintf(stderr, "expose event fwh=%d,%d, cwh=%d,%d\n", f->w, f->h, e->c->w, e->c->h);

	// workspace: dark blue background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 100 * (2 == i%3);

	struct camera_t *c = e->c;
	camera_capture(c);

	uint8_t *rgb = yuyv2rgb(c->head.start, c->w, c->h);
	//fillrgb(f->rgb

	int ox = (f->w - c->w) / 2;
	int oy = (f->h - c->h) / 2;
	assert(ox >= 0);
	assert(oy >= 0);
	for (int j = 0; j < (int)c->h; j++)
	for (int i = 0; i < (int)c->w; i++)
	for (int k = 0; k < 3; k++)
		f->rgb[3*(f->w*(j+oy)+(i+ox))+k] = rgb[3*(j*c->w+i)+k];
		//f->rgb[3*(f->w*(j+oy)+(i+ox))+k] = rgb[3*(j*c->w+(c->w -i-1))+k];
	free(rgb);
	f->changed = 1;
}

void kam_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);
	fprintf(stderr, "KEY k=%d ('%c') m=%d x=%d y=%d\n", k, k, m, x, y);

	// if ESC or q, exit
	if  (k == '\033' || k == 'q')
		ftr_notify_the_desire_to_stop_this_loop(f, 1);

	if (k == 'e')
		f->changed = 1;
}

int main()
{
	// state of the program
	struct kam_state e[1];

	// camera stuff
	int w = 800;
	int h = 600;
	e->c = camera_open("/dev/video0", w, h); // XXX: depends on cam!
	//e->c = camera_open("/dev/video0", 800, 448); // XXX: depends on cam!
	camera_init(e->c);
	camera_start(e->c);


	// window stuff
	struct FTR f = ftr_new_window(w+100,h+100);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , kam_key_handler);
	ftr_set_handler(&f, "idle", kam_exposer);
	int r = ftr_loop_run(&f);
	ftr_close(&f);

	camera_stop(e->c);
	camera_finish(e->c);
	camera_close(e->c);

	return 0;
}
