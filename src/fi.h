//
// FANCY IMAGE //
// -----------
//
// A data structure for dealing with large images

// The data structure.
//
// The fieds of this struct are intended to be read by the users.
struct fancy_image {
	int w;  // width
	int h;  // height
	int pd; // pixel dimension (samples per pixel)
	int no; // number of octaves
	char pad[200];
};


///////////////
// BASIC API //
///////////////

// open an image with the desired amount of cache
// (the cache size is honored only for tiled tiffs)
struct fancy_image fancy_image_open(const char *filename, double megabytes);

// close an image and free its associated cache
void fancy_image_close(struct fancy_image *f);

// obtain a sample of the image at the given point
// (if the point is outside the image domain, return NAN)
float fancy_image_getsample(struct fancy_image *f, int i,int j, int l);

// obtain a sample of the image at the given point and octave
// (if the coordinates or the octave are out of range, return NAN)
float fancy_image_getsample_oct(struct fancy_image *f,
		int octave, int i,int j, int l);




// UTILITY API: functions that have a more comfortable interface
// All these utility functions are defined in terms of the Basic Api above


// fill-in the array "out" with "f->pd" numbers
void fancy_image_getpixel(float *out, struct fancy_image *f, int i,int j);


// fill-in the array "out" with "f->pd" numbers
void fancy_image_getpixel_oct(float *out, struct fancy_image *f, int octave,
		int i,int j);

// getpixel with automatic scale selectoin, and trilinear interpolation
// and automatic choice of the octave
void fancy_image_trilinear(float *out, struct fancy_image *f,
		double x, double y, double dx, double dy);
