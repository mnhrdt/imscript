// TILOCT: a C library for reading small and large images
//
// A small image is an image that fits confortably in RAM memory.
// A large image is an image too large to fit in memory.
//
// This library provides an interface for reading transparently
// small and large images in a variety of formats.  Small images
// can be read from any file format (PNG, JPEG, TIFF, PNM...).
// Large images can only be read from tiled TIFF files.
//
//
// // BASIC API
// struct tiloct { int w, h, noctaves, spp, bps, fmt; ... };
// void   tiloct_read_image (struct tiloct *t, char *filename, int megabytes);
// void  *tiloct_get_pixel  (struct tiloct *t, int i, int j, int o);
//
// // CONVENIENCE API
// void  *tiloct_get_tile   (struct tiloct *t, int i, int j, int *ox, int *oy);
// void   tiloct_mark_read  (struct tiloct *t, int i, int j);
// void   tiloct_save_image (struct tiloct *t, char *filename);
//
// // FLOAT ACCESS
// void tiloct_get_pixel_float(float *out, struct tiloct *t, int i,int j,int o);
// float *tiloct_get_tile_float(struct tiloct *t, int i,int j, int *ox,int *oy);
//

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>


// 2^25 zoom factors should be enough for everybody
#define MAX_OCTAVES 25

struct tiloct {
	// user-visible data
	int w, h;     // image width, height
	int spp;      // samples per pixel
	int bps;      // bits per sample
	int fmt;      // sample format
	int noctaves; // number of octaves

	// tile geometry info
	int ntiles;   // total number of tiles (0 if not tiled)
	int tw, th;   // tile width, height
	int ta, td;   // tiles across, down

	char filename[MAX_OCTAVES][FILENAME_MAX];
};

//// visible interface to the library
//struct tiloct {
//	// data
//	int w, h;   // image width, height
//	int spp;    // samples per pixel
//	int bps;    // bits per sample
//	int fmt;    // sample format
//	int ntiles; // total number of tiles (0 if not tiled)
//
//	// only useful when tiled
//	int tw, th; // tile width, height
//	int ta, td; // tiles across, down
//};
//
//// a cache of tiles across several octaves
//struct tiff_octaves {
//	// user-visible data
//	int w, h, noctaves;
//	char filename[MAX_OCTAVES][FILENAME_MAX];
//	struct tiff_info i[MAX_OCTAVES];
//	void **c[MAX_OCTAVES];        // pointers to cached tiles
//
//	// data only necessary for garbage collection
//	unsigned int *a[MAX_OCTAVES]; // access counter for each tile
//	unsigned int ax; // global access counter
//	int curtiles;    // current number of tiles in memory
//	int maxtiles;    // number of tiles allowed to be in memory at once
//};


// how many units "u" are needed to cover a length "n"
static int how_many(int n, int u)
{
	return n/u + (_Bool)(n%u);

}

// read information from an open tiff file
static void get_tiff_info(struct tiff_info *t, TIFF *tif)
{
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &t->w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &t->h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &t->spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &t->bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &t->fmt);

	if (TIFFIsTiled(tif)) {
		TIFFGetField(tif, TIFFTAG_TILEWIDTH,  &t->tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &t->th);
		t->ta = how_many(t->w, t->tw);
		t->td = how_many(t->h, t->th);
		t->ntiles = TIFFNumberOfTiles(tif);
	}
}

// read information from a named tiff file
static int get_tiff_info_filename(struct tiff_info *t, char *fname)
{
	TIFF *tif = TIFFOpen(fname, "r");
	if (!tif) return 0;
	get_tiff_info(t, tif);
	TIFFClose(tif);
	return 1;
}

// read tile data from filename and pixel position
static void *read_tile_from_file(char *fname, int i, int j)
{
	TIFF *tif = TIFFOpen(fname, "r");
	struct tiff_info tinfo[1];
	get_tiff_info(tinfo, tif);
	void *r = malloc(TIFFTileSize(tif));
	TIFFReadTile(tif, r, i, j, 0, 0);
	TIFFClose(tif);
	return r;
}

// initialize a tile cache
void tiff_octaves_init(struct tiff_octaves *t, char *filepattern, int megabytes)
{
	// create filenames until possible
	t->noctaves = 0;
	for (int o = 0; o < MAX_OCTAVES; o++) {
		snprintf(t->filename[o], FILENAME_MAX, filepattern, o);
		if (!get_tiff_info_filename(t->i + o, t->filename[o])) break;
		if (o > 0 && 0 == strcmp(t->filename[o], t->filename[0])) break;
		t->noctaves += 1;

	}

	// set up essential data
	for (int o = 0; o < t->noctaves; o++) {
		t->c[o] = malloc(t->i[o].ntiles * sizeof*t->c);
		for (int j = 0; j < t->i[o].ntiles; j++)
			t->c[o][j] = 0;
	}

	// set up data for old tile deletion
	if (megabytes) {
		for (int o = 0; o < t->noctaves; o++)
			t->a[o] = malloc(t->i[o].ntiles * sizeof*t->a[o]);
		t->ax = 0;
		int tilesize = t->i->tw * t->i->th * (t->i->bps/8) * t->i->spp;
		t->maxtiles = megabytes * 1024.0 * 1024.0 / tilesize;
		t->curtiles = 0;
	} else
		t->a[0] = NULL;
}

// free a tile cache
void tiff_octaves_free(struct tiff_octaves *t)
{
	for (int i = 0; i < t->noctaves; i++) {
		for (int j = 0; j < t->i[i].ntiles; j++)
			free(t->c[i][j]);
		free(t->c[i]);
	}
}

// find the least recently accessed tile of a tile cache
static void find_oldest_tile(struct tiff_octaves *t, int *out_oct, int *out_idx)
{
	int omin = -1, imin = -1;
	for (int o = 0; o < t->noctaves; o++)
		for (int i = 0; i < t->i[o].ntiles; i++)
			if (t->a[o][i]) {
				if (imin >= 0) {
					if (t->a[o][i] < t->a[omin][imin]) {
						imin = i;
						omin = o;
					}
				} else {
					imin = i;
					omin = o;
				}
			}
	*out_idx = imin;
	*out_oct = omin;
}

// free the oldest (least recently accessed) tile of a tile cache
static void free_oldest_tile(struct tiff_octaves *t)
{
	int o, i;
	find_oldest_tile(t, &o, &i);

	free(t->c[o][i]);
	t->c[o][i] = 0;
	t->a[o][i] = 0;
	t->curtiles -= 1;
}

// notify that a tile has been accessed
static void notify_tile_access(struct tiff_octaves *t, int o, int i)
{
	t->a[o][i] = ++t->ax;
}

// enforce the condition  a <= x <= b
static int bound(int a, int x, int b)
{
	if (x < a) x = a;
	if (x > b) x = b;
	return x;
}

// compute the tile index of the given pixel
static int which_tile(struct tiff_info *t, int i, int j)
{
	if (i < 0 || j < 0 || i >= t->w || j >= t->h)
		return -1;
	int ti = i / t->tw;
	int tj = j / t->th;
	return tj * t->ta + ti;
}

// pointer to the tile that contains pixel (i,j) from octave o
void *tiff_octaves_gettile(struct tiff_octaves *t, int o, int i, int j)
{
	// sanitize input
	o = bound(0, o, t->noctaves - 1);
	i = bound(0, i, t->i[o].w - 1);
	j = bound(0, j, t->i[o].h - 1);

	// get valid tile index
	int tidx = which_tile(t->i + o, i, j);
	if (tidx < 0) return NULL;

	// if tile data does not yet exist, read it from file
	if (!t->c[o][tidx]) {
		if (t->a[0] && t->curtiles == t->maxtiles)
			free_oldest_tile(t);
		t->c[o][tidx] = read_tile_from_file(t->filename[o], i, j);
		t->curtiles += 1;
	}
	if (t->a[0])
		notify_tile_access(t, o, tidx);

	return t->c[o][tidx];
}

// pointer to the pixel (i,j) from octave o
void *tiff_octaves_getpixel(struct tiff_octaves *t, int o, int i, int j)
{
	void *tile = tiff_octaves_gettile(t, o, i, j);
	if (!tile) return NULL;

	// get pointer to requested pixel
	struct tiff_info *ti = t->i + o;
	int ii = i % ti->tw;
	int jj = j % ti->th;
	int pixel_index = jj * ti->tw + ii;
	int pixel_position = pixel_index * ti->spp * (ti->bps / 8);
	return pixel_position + (char*)tile;
}


// silly example
static int main_octaves(int c, char *v[])
{
	// read input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s inpattern npixels\n", *v);
		//                          0 1         2
		return 1;
	}
	char *filepattern = v[1];
	int npixels = atoi(v[2]);

	// initialize state
	int megabytes = 100;
	struct tiff_octaves t[1];
	tiff_octaves_init(t, filepattern, megabytes);

	// access some random pixels
	for (int i = 0; i < npixels; i++) {
		int o = rand() % t->noctaves;
		int p = rand() % t->i[o].w;
		int q = rand() % t->i[o].h;
		void *pixel = tiff_octaves_getpixel(t, o, p, q);
		fprintf(stderr, "p_%-3d  (%7d , %7d) = %p\n", o, p, q, pixel);
	}

	// cleanup and exit
	tiff_octaves_free(t);
	return 0;
}


#ifdef MAIN_OCTAVES
int main(int c, char *v[])
{
	return main_octaves(c, v);
}
#endif
