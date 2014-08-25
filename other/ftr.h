#ifndef _FTR_H
#define _FTR_H

// data structure to store the state of a window
struct FTR {
	// visible state
	int w, h;
	unsigned char *rgb;
	int changed;

	void *userdata; // ignored by the library

	// hidden implementation details
	char pad[200];
};

// type of a handler function
typedef void (*ftr_event_handler_t)(struct FTR*,int,int,int,int);

// window management
struct FTR ftr_new_window(int w, int h);
struct FTR ftr_new_window_with_image_uint8_rgb(unsigned char *i, int w, int h);
void ftr_fork_window_with_image_uint8_rgb(unsigned char *i, int w, int h);
void ftr_change_title(struct FTR *f, char *title);
void ftr_close(struct FTR *f);

// blocking calls
void ftr_wait_for_mouse_click(struct FTR *f, int *x, int *y);
void ftr_wait_for_mouse_click3(struct FTR *f, int *x, int *y, int *b);
void ftr_wait_for_key_depress(struct FTR *f, int *x, int *y, int *k, int *m);

// default and example handlers
void ftr_handler_exit_on_ESC(struct FTR*,int,int,int,int);
void ftr_handler_exit_on_ESC_or_q(struct FTR*,int,int,int,int);
void ftr_handler_toggle_idle(struct FTR*,int,int,int,int);
void ftr_handler_stop_loop(struct FTR*,int,int,int,int);
void ftr_handler_dummy(struct FTR*,int,int,int,int);

// event loop
int ftr_loop_run(struct FTR *f); // returns when the loop is finished
void ftr_loop_fork(struct FTR *f); // returns immediately, forks a new process
int ftr_num_pending(struct FTR *f);
void ftr_notify_the_desire_to_stop_this_loop(struct FTR *f, int return_value);

int ftr_set_handler(struct FTR *f, char *id, ftr_event_handler_t e);
ftr_event_handler_t ftr_get_handler(struct FTR *f, char *id);

// forking interface (so far not implemented)
//void ftr_signal_quit(struct FTR *);
//void ftr_signal_update(struct FTR *);


// ascii keys with name (necessary because '\e' is not standard)
#define FTR_KEY_ESC        27
#define FTR_KEY_DEL        127

// non-ascii keys (numbers inspired by glut)
#define FTR_KEY_FN         1000
#define FTR_KEY_LEFT       1100
#define FTR_KEY_UP         1101
#define FTR_KEY_RIGHT      1102
#define FTR_KEY_DOWN       1103
#define FTR_KEY_PAGE_UP    1104
#define FTR_KEY_PAGE_DOWN  1105
#define FTR_KEY_HOME       1106
#define FTR_KEY_END        1107
#define FTR_KEY_INSERT     1108

// key modifiers (numbers inspired by X)
#define FTR_MASK_SHIFT     1
#define FTR_MASK_LOCK      2
#define FTR_MASK_CONTROL   4
#define FTR_MASK_MOD1      8
#define FTR_MASK_MOD2      16
#define FTR_MASK_MOD3      32
#define FTR_MASK_MOD4      64
#define FTR_MASK_MOD5      128

// button modifiers (inspired by X)
#define FTR_BUTTON_LEFT    256
#define FTR_BUTTON_MIDDLE  512
#define FTR_BUTTON_RIGHT   1024
#define FTR_BUTTON_UP      2048
#define FTR_BUTTON_DOWN    4096

#endif//_FTR_H
