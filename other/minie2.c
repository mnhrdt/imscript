#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	ftr_set_handler(&f, "idle", ftr_handler_stop_loop);
	ftr_loop_run(&f);
	return printf("finí!\n");
}
