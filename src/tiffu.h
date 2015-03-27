#include <stdint.h>
#include <stdbool.h>

#ifndef FILENAME_MAX
#define FILENAME_MAX 1000
#endif

struct tiff_info {
	int32_t w;   // image width
	int32_t h;   // image height
	int16_t spp; // samples per pixel
	int16_t bps; // bits per sample
	int16_t fmt; // sample format
	bool broken; // whether pixels are contiguous or broken
	bool packed; // whether bps=1,2 or 4
	bool tiled;  // whether data is organized into tiles
	bool compressed;
	int ntiles;

	// only if tiled
	int32_t tw;  // tile width
	int32_t th;  // tile height
	int32_t ta;  // tiles across
	int32_t td;  // tiles down
};

struct tiff_tile_cache {
	// essential data
	//
	char filename[FILENAME_MAX];
	struct tiff_info i[1];
	void **c;        // pointers to cached tiles

	// data only necessary to delete tiles when the memory is full
	//
	unsigned int *a; // access counter for each tile
	unsigned int ax; // global access counter
	int curtiles;    // current number of tiles in memory
	int maxtiles;    // number of tiles allowed to be in memory at once
};
void tiff_tile_cache_init(struct tiff_tile_cache *t, char *fname, int mbytes);
void tiff_tile_cache_free(struct tiff_tile_cache *t);
void *tiff_tile_cache_getpixel(struct tiff_tile_cache *t, int i, int j);
void convert_pixel_to_float(float *out, struct tiff_info *t, void *in);
