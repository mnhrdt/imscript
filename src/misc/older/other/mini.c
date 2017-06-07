#include "ftr.c"

int main(void)
{
	struct FTR f = ftr_new_window(320, 200);
	return ftr_loop_run(&f);
}
