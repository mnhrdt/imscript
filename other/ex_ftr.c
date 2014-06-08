// port of cimg tutorial example to ftr

#include "ftr.h"

int main() {
	int w, h;
	unsigned char *x = ex_read_image_rgb24("/tmp/lenak.png", &w, &h);
	ex_inline_blur(x, w, h, 2.5);

	// OBS: "fork window" creates a new process with shared memory.
	// In unix and os/x: mmap(shared, anonymous) + fork()
	// In windows: CreateFileMapping and RunProcess
	struct FTR fx = ftr_fork_window_with_image(x, w, h);
	struct FTR fg = ftr_fork_window_with_image(0, 500, 400);

	while(ftr_is_running(&fx) && ftr_is_running(&fy))
	{
		ftr_wait(&fx);
		int b, y;
		if (ftr_had_button(&fx, &b, 0, 0, &y) && y >= 0 && y < h)
		{
			unsigned char *line = x + 3*w*y;
			ex_draw_rgb_graphs(fg.rgb, fg.w, fg.h, line, w);
			ftr_notify_change(&fg);
		}

	}
}
