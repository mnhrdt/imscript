struct FTR ftr_new_window(int w, int h)
{
	return ftr_new_window_with_image_uint8_rgb(NULL, w, h);
}

void ftr_fork_window_with_image_uint8_rgb(unsigned char *x, int w, int h)
{
	if (!fork())
	{
		struct FTR f = ftr_new_window_with_image_uint8_rgb(x, w, h);
		exit(ftr_loop_run(&f));
	}
}


void ftr_loop_fork(struct FTR *f)
{
	if (!fork())
		exit(ftr_loop_run(f));
}

void ftr_handler_exit_on_ESC(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_exit_on_ESC_or_q(struct FTR *f, int k, int m, int x, int y)
{
	if  (k == '\033' || k=='q' || k=='Q')
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_stop_loop(struct FTR *f, int k, int m, int x, int y)
{
	ftr_notify_the_desire_to_stop_this_loop(f, 0);
}

void ftr_handler_dummy(struct FTR *f, int k, int m, int x, int y)
{
}

void ftr_handler_toggle_idle(struct FTR *ff, int k, int m, int x, int y)
{
	struct _FTR *f = (void*)ff;
	ftr_event_handler_t tmp = f->handle_idle;
	f->handle_idle = f->handle_idle_toggled;
	f->handle_idle_toggled = tmp;
}

int ftr_set_handler(struct FTR *ff, char *id, ftr_event_handler_t e)
{
	struct _FTR *f = (void*)ff;
	if (0) ;
	else if (0 == strcmp(id, "key"))    { f->handle_key    = e; return 0; }
	else if (0 == strcmp(id, "button")) { f->handle_button = e; return 0; }
	else if (0 == strcmp(id, "motion")) { f->handle_motion = e; return 0; }
	else if (0 == strcmp(id, "expose")) { f->handle_expose = e; return 0; }
	else if (0 == strcmp(id, "resize")) { f->handle_resize = e; return 0; }
	else if (0 == strcmp(id, "idle"))   { f->handle_idle   = e; return 0; }
	return fprintf(stderr, "WARNING: unrecognized event \"%s\"\n", id);
}

ftr_event_handler_t ftr_get_handler(struct FTR *ff, char *id)
{
	struct _FTR *f = (void*)ff;
	if (0) ;
	else if (0 == strcmp(id, "key"))    return f->handle_key   ;
	else if (0 == strcmp(id, "button")) return f->handle_button;
	else if (0 == strcmp(id, "motion")) return f->handle_motion;
	else if (0 == strcmp(id, "expose")) return f->handle_expose;
	else if (0 == strcmp(id, "resize")) return f->handle_resize;
	else if (0 == strcmp(id, "idle"))   return f->handle_idle  ;
	return fprintf(stderr, "WARNING: bad event \"%s\"\n", id),
		(ftr_event_handler_t)NULL;
}

static void handle_click_wait(struct FTR *f, int b, int m, int x, int y)
{
	if (b == FTR_BUTTON_LEFT || b == FTR_BUTTON_RIGHT)
		ftr_notify_the_desire_to_stop_this_loop(f, 10000*y + x);
}

static void handle_click_wait3(struct FTR *f, int b, int m, int x, int y)
{
	if (b == FTR_BUTTON_LEFT || b == FTR_BUTTON_RIGHT)
	{
		int bit = b == FTR_BUTTON_LEFT;
		ftr_notify_the_desire_to_stop_this_loop(f, 2*(10000*y + x)+bit);
	}
}

void ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y)
{
	ftr_set_handler(f, "button", handle_click_wait);
	int r = ftr_loop_run(f);
	if (x) *x = r % 10000;
	if (y) *y = r / 10000;
}

void ftr_wait_for_mouse_click3(struct FTR *f, int *x, int *y, int *b)
{
	ftr_set_handler(f, "button", handle_click_wait3);
	int r = ftr_loop_run(f);
	int bit = r % 2;
	r /= 2;
	if (x) *x = r % 10000;
	if (y) *y = r / 10000;
	if (b) *b = bit ? FTR_BUTTON_LEFT : FTR_BUTTON_RIGHT;
}
