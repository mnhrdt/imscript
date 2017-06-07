#include "ftr.c"

void print_event(struct FTR *f, int k, int m, int x, int y)
{
	printf("event k=%d m=%d x=%d y=%d\n", k, m, x, y);
}

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	ftr_set_handler(&f, "key", print_event);
	ftr_set_handler(&f, "button", print_event);
	ftr_set_handler(&f, "motion", print_event);
	ftr_set_handler(&f, "resize", print_event);
	return ftr_loop_run(&f);
}
