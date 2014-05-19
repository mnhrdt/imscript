// data structure to store the state of a window
struct FTR;
// members: .w, .h, .rgb8_buffer

// type of a handler function
typedef void (*ftr_event_handler_t)(struct FTR*,int,int,int,int);


// window management
struct FTR ftr_new_window(void);
struct FTR ftr_new_window_and_malloc_buffer(int w, int w);
struct FTR ftr_new_window_from_rgba32_image(unsigned char *buf, int w, int w);
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
