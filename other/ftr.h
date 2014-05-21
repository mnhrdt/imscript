// data structure to store the state of a window
struct FTR;
// members: .w, .h, .rgb8_buffer

// type of a handler function
typedef void (*ftr_event_handler_t)(struct FTR*,int,int,int,int);


// window management
struct FTR ftr_new_window(void);
struct FTR ftr_new_window_and_malloc_buffer(int w, int h);
struct FTR ftr_new_window_with_image_uint8_argb(unsigned char *i, int w, int h);
void ftr_change_title(struct FTR *f, char *title);
void ftr_close(struct FTR *f);

// blocking calls
void ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y, int *b, int *m);
void ftr_wait_for_key_depress(struct FTR *f, int *x, int *y, int *k, int *m);

// event loop
int ftr_loop_run(struct FTR *f); // returns when the loop is finished
int ftr_set_handler(struct FTR *f, char *id, ftr_event_handler_t e);
//int ftr_loop_fork(struct FTR *f); // returns immediately

//// default and example handlers
//ftr_event_handler ftr_handler_exit_on_ESC;

// convenience functions
//void ftr_fork_window_with_rgba32_image(unsigned char *buf, int w, int h);
//void ftr_fork_window_with_rgb32_image(unsigned char *buf, int w, int h);
//void ftr_fork_window_with_f32_image(float *buf, int w, int h);
//void ftr_fork_window_with_f32_image_vec(float *buf, int w, int h, int pd);

struct FTR {
	// state
	int w, h, max_w, max_h;
	unsigned char *argb;
	int do_exit;
	int changed;

	// user-supplied handlers
	ftr_event_handler_t handle_key;
	ftr_event_handler_t handle_button;
	ftr_event_handler_t handle_motion;
	ftr_event_handler_t handle_expose;
	ftr_event_handler_t handle_resize;
	ftr_event_handler_t handle_idle;


	char pad[100];
	//// X11-related internal data
	//Display *display;
	//Visual *visual;
	//Window window;
	//GC gc;
	//XImage *ximage;
	//int imgupdate;
};
