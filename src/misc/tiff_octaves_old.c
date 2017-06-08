#include <assert.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <tiffio.h>



// 2^25 zoom factors should be enough for everybody
#define MAX_OCTAVES 25




// tiff file information (without the actual data)
struct tiff_info {
	int32_t w;       // image width
	int32_t h;       // image height
	int16_t spp;     // samples per pixel
	int16_t bps;     // bits per sample
	int16_t fmt;     // sample format
	bool broken;     // whether pixels are contiguous or broken
	bool packed;     // whether bps=1,2 or 4
	bool tiled;      // whether data is organized into tiles
	bool compressed; // whether the tile data is compressed
	int ntiles;      // total number of tiles (0 if not tiled)

	// only if tiled
	int32_t tw;      // tile width
	int32_t th;      // tile height
	int32_t ta;      // tiles across
	int32_t td;      // tiles down
};

// tiff tile with a pointer to the data
struct tiff_tile {
	int32_t w, h;
	int16_t spp, bps, fmt;
	uint8_t *data;
};

// a cache of tiles across several octaves
struct tiff_octaves {
	// essential data
	//
	int noctaves;
	char filename[MAX_OCTAVES][FILENAME_MAX];
	struct tiff_info i[MAX_OCTAVES];
	void **c[MAX_OCTAVES];        // pointers to cached tiles

	// data only necessary for garbage collection
	//
	unsigned int *a[MAX_OCTAVES]; // access counter for each tile
	unsigned int ax; // global access counter
	int curtiles;    // current number of tiles in memory
	int maxtiles;    // number of tiles allowed to be in memory at once
};



// print an error message and abort the program
static void fail(const char *fmt, ...)
{
	va_list argp;
	fprintf(stderr, "\nFAIL: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	abort();
}

// like malloc, but always return a valid pointer
static void *xmalloc(size_t size)
{
	if (size == 0)
		fail("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
		fail("xmalloc: out of memory requesting %zu bytes", size);
	return new;
}

// how many units "u" are needed to cover a length "n"
static int how_many(int n, int u)
{
	return n/u + (bool)(n%u);

}

// read information from an open tiff file
static void get_tiff_info(struct tiff_info *t, TIFF *tif)
{
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &t->w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &t->h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &t->spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &t->bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &t->fmt);

	uint16_t planarity, compression;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG,    &planarity);
	TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION,     &compression);
	t->broken = planarity != PLANARCONFIG_CONTIG;
	t->compressed = compression != COMPRESSION_NONE;
	t->packed = 0 != t->bps % 8;
	t->tiled = TIFFIsTiled(tif);

	if (t->tiled) {
		TIFFGetField(tif, TIFFTAG_TILEWIDTH,  &t->tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &t->th);
		t->ta = how_many(t->w, t->tw);
		t->td = how_many(t->h, t->th);
		t->ntiles = TIFFNumberOfTiles(tif);
		assert(t->ta * t->td == t->ntiles);
	}
}

// read information from a named tiff file
static bool get_tiff_info_filename(struct tiff_info *t, const char *fname)
{
	TIFF *tif = TIFFOpen(fname, "r");
	if (!tif)
		return false;
	get_tiff_info(t, tif);
	TIFFClose(tif);
	return true;
}

// read tile data from filename and pixel position
static void read_tile_from_file(struct tiff_tile *t, char *fname, int i, int j)
{
	TIFF *tif = TIFFOpen(fname, "r");
	if (!tif) fail("could not open TIFF file \"%s\" for reading", fname);

	struct tiff_info tinfo[1];
	get_tiff_info(tinfo, tif);

	if (tinfo->broken) fail("broken pixels not supported yet");
	if (!tinfo->tiled) fail("I do not read non-tiled tiffs");

	int tbytes = TIFFTileSize(tif);
	t->data = xmalloc(tbytes);
	memset(t->data, 0, tbytes);
	tsize_t r = TIFFReadTile(tif, t->data, i, j, 0, 0);
	if (r == -1) fail("could not read a tile from file \"%s\"", fname);

	TIFFClose(tif);
}

// initialize a tile cache
void tiff_octaves_init(struct tiff_octaves *t, const char *filepattern, int megabytes)
{
	// create filenames until possible
	t->noctaves = 0;
	for (int o = 0; o < MAX_OCTAVES; o++)
	{
		snprintf(t->filename[o], FILENAME_MAX, filepattern, o);
		if (!get_tiff_info_filename(t->i + o, t->filename[o])) break;
		if (t->i[o].bps < 8 || t->i[o].packed)
			fail("caching of packed samples is not supported");
		if (o > 0) { // check consistency
			if (0 == strcmp(t->filename[o], t->filename[0])) break;
			if (t->i[o].bps != t->i->bps) fail("inconsistent bps");
			if (t->i[o].spp != t->i->spp) fail("inconsistent spp");
			if (t->i[o].fmt != t->i->fmt) fail("inconsistent fmt");
			if (t->i[o].tw != t->i->tw) fail("inconsistent tw");
			if (t->i[o].th != t->i->th) fail("inconsistent th");
		}
		t->noctaves += 1;

	}
	if (t->noctaves < 1)
		fail("could not get any file with pattern \"%s\"", filepattern);

	// set up essential data
	for (int o = 0; o < t->noctaves; o++)
	{
		t->c[o] = xmalloc(t->i[o].ntiles * sizeof*t->c);
		for (int j = 0; j < t->i[o].ntiles; j++)
			t->c[o][j] = 0;
	}

	// set up data for old tile deletion
	if (megabytes) {
		for (int o = 0; o < t->noctaves; o++)
			t->a[o] = xmalloc(t->i[o].ntiles * sizeof*t->a[o]);
		t->ax = 0;
		int tilesize = t->i->tw * t->i->th * (t->i->bps/8) * t->i->spp;
		double mbts = tilesize / (1024.0 * 1024);
		t->maxtiles = megabytes / mbts;
		t->curtiles = 0;
	} else  {
		// unlimited tile usage
		t->a[0] = NULL;
	}
}

// free a tile cache
void tiff_octaves_free(struct tiff_octaves *t)
{
	for (int i = 0; i < t->noctaves; i++)
	{
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
	assert(imin >= 0);
	assert(t->a[omin][imin] > 0);

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
	int r = tj * t->ta + ti;
	if (r < 0 || r >= t->ntiles)
		fail("bad tile index %d for point (%d %d)", r, i, j);
	return r;
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
		struct tiff_tile tmp[1];
		read_tile_from_file(tmp, t->filename[o], i, j);
		t->c[o][tidx] = tmp->data;
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

// convenience function
void tiff_octaves_float_tile(float *out, struct tiff_octaves *t,
							int o, int i, int j)
{
	void *in = tiff_octaves_gettile(t, o, i, j);

	struct tiff_info *ti = t->i + o;
	int tsamples = ti->tw * ti->th * ti->spp;
	for (int l = 0; l < tsamples; l++) {
		switch(ti->fmt) {
		case SAMPLEFORMAT_UINT:
			if (8 == ti->bps)        out[l] = ((uint8_t *)in)[l];
			else if (16 == ti->bps)  out[l] = ((uint16_t *)in)[l];
			else if (32 == ti->bps)  out[l] = ((uint32_t *)in)[l];
			break;
		case SAMPLEFORMAT_INT:
			if (8 == ti->bps)        out[l] = ((int8_t *)in)[l];
			else if (16 == ti->bps)  out[l] = ((int16_t *)in)[l];
			else if (32 == ti->bps)  out[l] = ((int32_t *)in)[l];
			break;
		case SAMPLEFORMAT_IEEEFP:
			if (32 == ti->bps)       out[l] = ((float *)in)[l];
			else if (64 == ti->bps)  out[l] = ((double *)in)[l];
			break;
		}

	}
}


#ifdef MAIN_OCTAVES
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
	for (int i = 0; i < npixels; i++)
	{
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


int main(int c, char *v[])
{
	return main_octaves(c, v);
}
#endif
