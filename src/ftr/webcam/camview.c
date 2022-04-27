#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ftr.h"
#include "cam.h"

#include "seconds.c"

struct kam_state {
	struct camera *c;
	uint8_t *prev;
};

static void kam_exposer(struct FTR *f, int b, int m, int unused_x, int unused_y)
{
	static double tic = 0;
	double tac = seconds();
	fprintf(stderr, "tictac = %g ms (%g fps)\n",1000*(tac-tic),1/(tac-tic));
	tic = tac;

	(void)unused_x; (void)unused_y;
	struct kam_state *e = f->userdata;

	// workspace: dark blue background
	for (int i = 0; i < 3 * f->w * f->h; i++)
		f->rgb[i] = 100 * (2 == i%3);

	struct camera *c = e->c;
	camera_grab_rgb(c);

	// leave a symmetric margin
	int ox = (f->w - c->w) / 2;
	int oy = (f->h - c->h) / 2;
	assert(ox >= 0);
	assert(oy >= 0);

	//for (int j = 0; j < c->h; j++)
	//	memcpy(f->rgb+ 3*(f->w*(j+oy)+ox), c->rgb+ 3*(c->w*j), 3*c->w);

	for (int j = 0; j < c->h; j++)
	for (int i = 0; i < c->w; i++)
	for (int k = 0; k < 3   ; k++)
	{
		float v1 = c->rgb [3*(j*c->w+i)+k];
		float v2 = e->prev[3*(j*c->w+i)+k];
		float v = 3*(v2 - v1) + 127;
		if (v > 255) v = 255;
		if (v < 0) v = 0;
		f->rgb[3*(f->w*(j+oy)+i+ox)+k] = v;
	}

	memcpy(e->prev, c->rgb, 3 * c->w * c->h);
}

void kam_key_handler(struct FTR *f, int k, int m, int x, int y)
{
	if (m & FTR_MASK_SHIFT && islower(k)) k = toupper(k);
	fprintf(stderr, "KEY k=%d ('%c') m=%d x=%d y=%d\n", k, k, m, x, y);

	if  (k == '\033' || k == 'q') // ESC or q
		ftr_notify_the_desire_to_stop_this_loop(f, 1);
}

int main()
{
	// state of the program
	struct kam_state e[1];

	// camera stuff
	int w = 640;//800;
	int h = 360;//600;
	e->c = camera_begin("/dev/video2", w, h); // XXX: depends on cam!

	e->prev = malloc(2*3 * w * h);

	// window stuff
	struct FTR f = ftr_new_window(w + 100, h + 100);
	f.userdata = e;
	f.changed = 1;
	ftr_set_handler(&f, "key"   , kam_key_handler);
	ftr_set_handler(&f, "idle", kam_exposer);
	int r = ftr_loop_run(&f);
	ftr_close(&f);

	free(e->prev);
	camera_end(e->c);


	return 0;
}
