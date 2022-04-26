#include <stdint.h> // uint8_t, uint32_t
#include <stdlib.h>
#include <sys/time.h> // timeval

// used only inside the camera struct
struct camera_buffer_t {
	uint8_t* start;
	size_t length;
};

// internal struct, just for this program
struct camera_t {
	// intended for external use
	uint32_t w;
	uint32_t h;
	uint8_t *rgb;

	// implementation details
	int fd;
	size_t buffer_count;
	struct camera_buffer_t *buffers;
	struct camera_buffer_t head;
};

// open the device file and fill-in struct fields
// returns an allocated struct that needs to be freed afterwards
struct camera_t* camera_open(const char *device, uint32_t w, uint32_t h);

// check and setup capabilities via ioctls, memmap buffers
// (the actual webcam is still off)
void camera_init(struct camera_t* c);

// switch on the camera and start recording frames into the buffers
void camera_start(struct camera_t* c);

// switch off the camera
void camera_stop(struct camera_t* c);

// free buffers and stuff
void camera_finish(struct camera_t* c);

// close device and free struct memory
void camera_close(struct camera_t* c);

// copy one frame from the buffer ring to the head, at any time
int camera_capture(struct camera_t* c);

// wait at most "timeout" and capture the next frame
int camera_frame(struct camera_t* c, struct timeval timeout);

uint8_t* yuyv2rgb(uint8_t* yuyv, uint32_t width, uint32_t height);
void fillrgb(uint8_t *rgb, uint8_t *yuyv, uint32_t width, uint32_t height);
