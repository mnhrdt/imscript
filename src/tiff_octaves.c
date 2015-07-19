#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <tiffio.h>


// structs {{{1

struct tiff_tile {
	int32_t w, h;
	int16_t spp, bps, fmt;
	bool broken;
	uint8_t *data;

	int offx, offy;
};

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

// general utility functions {{{1


// exit the program printing an error message
#ifndef _FAIL_C
#define _FAIL_C
static void fail(const char *fmt, ...)
{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#  ifdef NDEBUG
	exit(-1);
#  else//NDEBUG
	exit(*(int *)0x43);
#  endif//NDEBUG
}
#endif//_FAIL_C


// like malloc, but always returns a valid pointer
#ifndef _XMALLOC_C
#define _XMALLOC_C
static void *xmalloc(size_t size)
{
	if (size == 0)
		fail("xmalloc: zero size");
	void *p = malloc(size);
	if (!p)
	{
		double sm = size / (0x100000 * 1.0);
		fail("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return p;
}
#endif//_XMALLOC_C


static int my_computetile(struct tiff_info *t, int i, int j)
{
	if (i < 0 || j < 0 || i >= t->w || j >= t->h)
		return -1;
		//fail("got bad pixel (%d %d) [%d %d]", i, j, t->w, t->h);
	int ti = i / t->tw;
	int tj = j / t->th;
	int r = tj * t->ta + ti;
	if (r < 0 || r >= t->ntiles)
		fail("bad tile index %d for point (%d %d)", r, i, j);
	return r;
}


// open a TIFF file, with some magic to access subimages
// (i.e., filename "file.tif,3" refers to the third sub-image)
static TIFF *tiffopen_fancy(char *filename, char *mode)
{
	//fprintf(stderr, "tiffopen fancy \"%s\",\"%s\"\n", filename, mode);
	char *comma = strrchr(filename, ',');
	if (*mode != 'r' || !comma)
	def:	return TIFFOpen(filename, mode);

	int aftercomma = strlen(comma + 1);
	int ndigits = strspn(comma + 1, "0123456789");

	if (aftercomma != ndigits) goto def;

	char buf[FILENAME_MAX];
	strncpy(buf, filename, FILENAME_MAX);
	comma = strrchr(buf, ',');
	*comma = '\0';
	int index = atoi(comma + 1);

	TIFF *tif = TIFFOpen(buf, mode);
	if (!tif) return tif;
	for (int i = 0; i < index; i++)
		TIFFReadDirectory(tif);
	
	return tif;
}


// tell how many units "u" are needed to cover a length "n"
static int how_many(int n, int u)
{
	assert(n > 0);
	assert(u > 0);
	return n/u + (bool)(n%u);

}



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
	} else {
		t->ta = t->td = 1;
		t->tw = t->w;
		t->th = t->h;
		t->ntiles = 1;
	}
}

static bool get_tiff_info_filename_e(struct tiff_info *t, char *fname)
{
	TIFF *tif = tiffopen_fancy(fname, "r");
	if (!tif)
		return false;
	get_tiff_info(t, tif);
	TIFFClose(tif);
	return true;
}

static int tiff_imagewidth(TIFF *tif)
{
	uint32_t w, r = TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	if (r != 1) fail("can not read TIFFTAG_IMAGEWIDTH");
	return w;
}

static int tiff_imagelength(TIFF *tif)
{
	uint32_t h, r = TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	if (r != 1) fail("can not read TIFFTAG_IMAGELENGTH");
	return h;
}

static int tiff_tilewidth(TIFF *tif)
{
	uint32_t R;
	int r = TIFFGetField(tif, TIFFTAG_TILEWIDTH, &R);
	if (r != 1) fail("TIFF tiles of unknown width");
	return R;
}

static int tiff_tilelength(TIFF *tif)
{
	uint32_t R;
	int r = TIFFGetField(tif, TIFFTAG_TILELENGTH, &R);
	if (r != 1) fail("TIFF tiles of unknown width");
	return R;
}

static int tiff_tilesacross(TIFF *tif)
{
	int w = tiff_imagewidth(tif);
	int tw = tiff_tilewidth(tif);
	return how_many(w, tw);
}

static int tiff_tilesdown(TIFF *tif)
{
	int h = tiff_imagelength(tif);
	int th = tiff_tilelength(tif);
	return how_many(h, th);
}

static int tiff_tile_corner(int p[2], TIFF *tif, int tidx)
{
	p[0] = p[1] = -1;

	int tw = tiff_tilewidth(tif);
	int th = tiff_tilelength(tif);
	int ta = tiff_tilesacross(tif);
	int td = tiff_tilesdown(tif);

	int i = tidx % ta;
	int j = tidx / ta;

	if (!(i < ta && j < td))
		return 0;

	p[0] = tw*i;
	p[1] = th*j;
	return 1;
}



static void read_scanlines(struct tiff_tile *tout, TIFF *tif)
{
	// read all file info
	struct tiff_info tinfo[1];
	get_tiff_info(tinfo, tif);

	// fill-in output information
	tout->w = tinfo->w;
	tout->h = tinfo->h;
	tout->bps = tinfo->bps;
	tout->fmt = tinfo->fmt;
	tout->spp = tinfo->spp;

	// define useful constants
	int pixel_size = tinfo->spp * tinfo->bps/8;
	int output_size = tout->w * tout->h * pixel_size;

	// allocate space for output data
	tout->data = xmalloc(output_size);

	// copy scanlines
	int scanline_size = TIFFScanlineSize(tif);
	assert(scanline_size == tinfo->w * pixel_size);
	for (int j = 0; j < tinfo->h; j++)
	{
		uint8_t *buf = tout->data + scanline_size*j;
		int r = TIFFReadScanline(tif, buf, j, 0);
		if (r < 0) fail("could not read scanline %d", j);
	}
}

static tsize_t my_readtile(TIFF *tif, tdata_t buf,
		uint32 x, uint32 y, uint32 z, tsample_t sample)
{
	tsize_t r = TIFFReadTile(tif, buf, x, y, z, sample);
	if (r == -1) memset(buf, 0, r = TIFFTileSize(tif));
	return r;
}

static void read_tile_from_file(struct tiff_tile *t, char *filename, int tidx)
{
	TIFF *tif = tiffopen_fancy(filename, "r");
	if (!tif) fail("could not open TIFF file \"%s\" for reading", filename);

	uint32_t w, h;
	uint16_t spp, bps, fmt, planarity;
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &fmt);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG,    &planarity);
	t->spp = spp;
	t->bps = bps;
	t->fmt = fmt;
	t->broken = false;

	if (planarity != PLANARCONFIG_CONTIG)
		fail("broken pixels not supported yet");

	if (TIFFIsTiled(tif)) {
		int nt = TIFFNumberOfTiles(tif);
		if (tidx < 0 || tidx >= nt)
			fail("bad tile index %d", tidx);

		t->w = tiff_tilewidth(tif);;
		t->h = tiff_tilelength(tif);

		int ii[2];
		tiff_tile_corner(ii, tif, tidx);
		//int ii = tw*i;
		//int jj = th*j;

		int tbytes = TIFFTileSize(tif);
		t->data = xmalloc(tbytes);
		memset(t->data, 0, tbytes);
		int r = my_readtile(tif, t->data, ii[0], ii[1], 0, 0);
		if (r != tbytes) fail("could not read tile");
	} else { // not tiled, read the whole image into 0th tile
		read_scanlines(t, tif);
	}

	TIFFClose(tif);
}




// getpixel cache with octaves {{{1

#define MAX_OCTAVES 25
struct tiff_octaves {
	// essential data
	//
	int noctaves;
	char filename[MAX_OCTAVES][FILENAME_MAX];
	struct tiff_info i[MAX_OCTAVES];
	void **c[MAX_OCTAVES];        // pointers to cached tiles

	// data only necessary to delete tiles when the memory is full
	//
	unsigned int *a[MAX_OCTAVES]; // access counter for each tile
	unsigned int ax; // global access counter
	int curtiles;    // current number of tiles in memory
	int maxtiles;    // number of tiles allowed to be in memory at once
};

//#include "smapa.h"
//SMART_PARAMETER(FIRST_OCTAVE,0)

void tiff_octaves_init(struct tiff_octaves *t, char *filepattern, int megabytes)
{
	fprintf(stderr, "tiff octaves init \"%s\"(%dMB)\n", filepattern, megabytes);
	// create filenames until possible
	t->noctaves = 0;
	for (int o = 0; o < MAX_OCTAVES; o++)
	{
		//int oo = o + FIRST_OCTAVE();
		snprintf(t->filename[o], FILENAME_MAX, filepattern, o);
		//fprintf(stderr, "f[%d]=%s\n", o, t->filename[o]);
		if (!get_tiff_info_filename_e(t->i + o, t->filename[o]))
			break;
		if (t->i[o].bps < 8 || t->i[o].packed)
			fail("caching of packed samples is not supported");
		if (0) {
			fprintf(stderr, "\tw = %d\n", (int)t->i[o].w);
			fprintf(stderr, "\th = %d\n", (int)t->i[o].h);
			fprintf(stderr, "\ttiled = %d\n", t->i[o].tiled);
			fprintf(stderr, "\ttw = %d\n", (int)t->i[o].tw);
			fprintf(stderr, "\tth = %d\n", (int)t->i[o].th);
		}
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
		fail("Could not get any file with pattern \"%s\"", filepattern);

	// set up essential data
	for (int o = 0; o < t->noctaves; o++)
	{
		t->c[o] = xmalloc((1 + t->i[o].ntiles) * sizeof*t->c);
		for (int j = 0; j < t->i[o].ntiles; j++)
			t->c[o][j] = 0;
	}

	// print debug info
	fprintf(stderr, "%d octaves:\n", t->noctaves);
	for (int o = 0; o < t->noctaves; o++)
	{
		struct tiff_info *ti = t->i + o;
		fprintf(stderr, "\toctave %d:", o);
		fprintf(stderr, " %dx%d", ti->w, ti->h);
		fprintf(stderr, " %d tiles (%dx%d) of size %dx%d",
				ti->ntiles, ti->ta, ti->td, ti->tw, ti->th);
		fprintf(stderr, "\n");
	}

	// set up data for old tile deletion
	if (megabytes) {
		for (int o = 0; o < t->noctaves; o++)
			t->a[o] = malloc(t->i[o].ntiles * sizeof*t->a[o]);
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

void tiff_octaves_free(struct tiff_octaves *t)
{
	for (int i = 0; i < t->noctaves; i++)
	{
		for (int j = 0; j < t->i[i].ntiles; j++)
			free(t->c[i][j]);
		free(t->c[i]);
	}
}

static int bound(int a, int x, int b)
{
	if (x < a) x = a;
	if (x > b) x = b;
	return x;
}

static void free_oldest_tile_octave(struct tiff_octaves *t)
{
	// find oldest tile
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

	// free it
	//
	//fprintf(stderr, "CACHE: FREEing tile %d of octave %d\n", imin, omin);
	free(t->c[omin][imin]);
	t->c[omin][imin] = 0;
	t->a[omin][imin] = 0;
	t->curtiles -= 1;
}

//static void free_oldest_half_of_tiles(struct tiff_octaves *t)
//{
//
//}

static void notify_tile_access_octave(struct tiff_octaves *t, int o, int i)
{
	t->a[o][i] = ++t->ax;
}

void *tiff_octaves_gettile(struct tiff_octaves *t, int o, int i, int j)
{
	// sanitize input
	o = bound(0, o, t->noctaves - 1);
	i = bound(0, i, t->i[o].w - 1);
	j = bound(0, j, t->i[o].h - 1);

	// get valid tile index
	int tidx = my_computetile(t->i + o, i, j);
	if (tidx < 0) return NULL;

	// if tile does not exist, read it from file
	if (!t->c[o][tidx])
//#pragma omp critical
	{
		if (t->a[0] && t->curtiles == t->maxtiles)
			free_oldest_tile_octave(t);

		fprintf(stderr,"CACHE: LOADing tile %d of octave %d\n",tidx,o);
		struct tiff_tile tmp[1];
		read_tile_from_file(tmp, t->filename[o], tidx);
		t->c[o][tidx] = tmp->data;

		t->curtiles += 1;
	}
	if (t->a[0])
		notify_tile_access_octave(t, o, tidx);

	return t->c[o][tidx];
}

void *tiff_octaves_getpixel(struct tiff_octaves *t, int o, int i, int j)
{
	//fprintf(stderr, "t_o_g(%d, %d, %d)\n", o, i, j);
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

// shy versions of the previous two functions
// (the shy functions avoid reading the disk)
void *tiff_octaves_gettile_shy(struct tiff_octaves *t, int o, int i, int j)
{
	// sanitize input
	o = bound(0, o, t->noctaves - 1);
	i = bound(0, i, t->i[o].w - 1);
	j = bound(0, j, t->i[o].h - 1);

	// get valid tile index
	int tidx = my_computetile(t->i + o, i, j);
	if (tidx < 0) return NULL;

	// if tile does not exist, return NULL
	if (!t->c[o][tidx]) return NULL;

	if (t->a[0])
		notify_tile_access_octave(t, o, tidx);

	return t->c[o][tidx];
}

void *tiff_octaves_getpixel_shy(struct tiff_octaves *t, int o, int i, int j)
{
	void *tile = tiff_octaves_gettile_shy(t, o, i, j);
	if (!tile) return NULL;

	// get pointer to requested pixel
	struct tiff_info *ti = t->i + o;
	int ii = i % ti->tw;
	int jj = j % ti->th;
	int pixel_index = jj * ti->tw + ii;
	int pixel_position = pixel_index * ti->spp * (ti->bps / 8);
	return pixel_position + (char*)tile;
}
