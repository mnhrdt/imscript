#include <stdint.h> // uint8_t
#include <stdlib.h>
#include <sys/time.h> // timeval

// internal struct, just for this program
struct camera {
	// intended for external use
	int w;
	int h;
	uint8_t *rgb;

	// implementation details
	int fd;
	size_t buffer_count;
	struct camera_buffer {
		uint8_t *start;
		size_t length;
	} *buffers;
	struct camera_buffer head;
};

// high-level API: just the public struct and the following three functions:
struct camera *camera_begin(const char *device, int w, int h);
void camera_end(struct camera *c);
void camera_grab_rgb(struct camera *c); // fill-in the c->rgb field


// open the device file and fill-in struct fields
// returns an allocated struct that needs to be freed afterwards
struct camera *camera_open(const char *device, int w, int h);

// check and setup capabilities via ioctls, memmap buffers
// (the actual webcam is still off)
void camera_init(struct camera *c);

// switch on the camera and start recording frames into the buffers
void camera_start(struct camera *c);

// switch off the camera
void camera_stop(struct camera *c);

// free buffers and stuff
void camera_finish(struct camera *c);

// close device and free struct memory
void camera_close(struct camera *c);

// copy one frame from the buffer ring to the head, at any time
int camera_capture(struct camera *c);

// wait at most "timeout" and capture the next frame
int camera_frame(struct camera *c, struct timeval timeout);

uint8_t *yuyv2rgb(uint8_t *yuyv, int width, int height);
void fillrgb(uint8_t *rgb, uint8_t *yuyv, int width, int height);
