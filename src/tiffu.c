// tiff utils
//
// info   f.tiff            # print misc info
// ntiles f.tiff            # print the total number of tiles
// tget   f.tiff n t.tiff   # get the nth tile (sizes must coincide)
// tput   f.tiff n t.tiff   # an image into the nth tile (sizes must coincide)
//
// meta   prog in.tiff ... out.tiff # run "prog" for all the tiles


#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>

#include "fail.c"
#include "xmalloc.c"


struct tiff_tile {
	int32_t w, h;
	int16_t spp, bps, fmt;
	bool broken;
	uint8_t *data;

	int offx, offy;
};

struct tiff_file_info {
	int32_t w;   // image width
	int32_t h;   // image height
	int16_t spp; // samples per pixel
	int16_t bps; // bits per sample
	int16_t fmt; // sample format
	bool broken; // whether pixels are contiguous or broken
	bool packed; // whether bps=1,2 or 4
	bool tiled;  // whether data is organized into tiles

	// only if tiled
	int32_t tw;  // tile width
	int32_t th;  // tile height
	int32_t ta;  // tiles across
	int32_t td;  // tiles down
};


// tell how many units "u" are needed to cover a length "n"
static int how_many(int n, int u)
{
	assert(n > 0);
	assert(u > 0);
	return n/u + (bool)(n%u);

}

static void get_tiff_file_info(struct tiff_file_info *t, TIFF *tif)
{
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &t->w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &t->h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &t->spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &t->bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &t->fmt);

	uint16_t planarity;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG,    &planarity);
	t->broken = planarity != PLANARCONFIG_CONTIG;
	t->packed = 0 != t->bps % 8;
	t->tiled = TIFFIsTiled(tif);

	if (t->tiled) {
		TIFFGetField(tif, TIFFTAG_TILEWIDTH,  &t->tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &t->th);
		t->ta = how_many(t->w, t->tw);
		t->td = how_many(t->h, t->th);
	}
}

static void get_tiff_file_info_filename(struct tiff_file_info *t, char *fname)
{
	TIFF *tif = TIFFOpen(fname, "r");
	if (!tif)
		fail("could not open TIFF file \"%s\" for reading", fname);
	get_tiff_file_info(t, tif);
	TIFFClose(tif);
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

static int tiff_samplesperpixel(TIFF *tif)
{
	uint16_t spp;
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	return spp;
}

// divide by 8 to obtain the size per sample (then, 0=packed data)
static int tiff_bitspersample(TIFF *tif)
{
	uint16_t bps;
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &bps);
	return bps;
}

static int tiff_sampleformat(TIFF *tif)
{
	uint16_t fmt;
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
	return fmt;
}

static bool tiff_brokenpixels(TIFF *tif)
{
	uint16_t planarity;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &planarity);
	return planarity != PLANARCONFIG_CONTIG;
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


// write a non-tiled TIFF image
static void write_tile_to_file(char *filename, struct tiff_tile *t)
{
	if (t->broken) fail("can not save broken tiles yet");

	TIFF *tif = TIFFOpen(filename, "w");
	if (!tif) fail("could not open TIFF file \"%s\" for writing", filename);

	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,      t->w);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH,     t->h);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, t->spp);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,   t->bps);
	TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,    t->fmt);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG,    PLANARCONFIG_CONTIG);

	int Bps = t->bps/8;
	if (!Bps) fail("packed pixels not supported yet");

	int scanline_size = t->w * t->spp * Bps;
	for (int i = 0; i < t->h; i++)
	{
		void *line = i*scanline_size + t->data;
		int r = TIFFWriteScanline(tif, line, i, 0);
		if (r < 0) fail("error writing %dth TIFF scanline", i);
	}

	TIFFClose(tif);
}

static void dalloc_tile_from_open_file(struct tiff_tile *t, TIFF *tif, int tidx)
{
	//struct tiff_file_info tinfo[1];
	//get_tiff_file_info(tinfo, tif);

	//int nt = TIFFNumberOfTiles(tif);
	//if (tidx < 0 || tidx >= nt)
	//	fail("bad tile index %d", tidx);

	//t->w = tinfo->tw;
	//t->h = tinfo->th;
	//t->spp = tinfo->spp;
	//t->bps = tinfo->bps;
	//t->fmt = tinfo->fmt;
	//t->broken = false;

	//int tbytes = TIFF

	//int i = tidx % ta;
	//int j = tidx / ta;

	//int ii = tw*i;
	//int jj = th*i;

	//int tbytes = TIFFTileSize(tif);
	//t->data = xmalloc(tbytes);
	//memset(t->data, 0, tbytes);
	//int r = TIFFReadTile(tif, t->data, ii, jj, 0, 0);
	//if (r != tbytes) fail("could not read tile");

}

static void read_tile_from_file(struct tiff_tile *t, char *filename, int tidx)
{
	TIFF *tif = TIFFOpen(filename, "r");
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
		fail("broken pixels not supporte yet");

	if (TIFFIsTiled(tif)) {
		int nt = TIFFNumberOfTiles(tif);
		if (tidx < 0 || tidx >= nt)
			fail("bad tile index %d", tidx);

		int tw = tiff_tilewidth(tif);
		int th = tiff_tilelength(tif);
		int ta = tiff_tilesacross(tif);
		int td = tiff_tilesdown(tif);

		t->w = tw;
		t->h = th;

		int i = tidx % ta;
		int j = tidx / ta;

		int ii = tw*i;
		int jj = th*j;

		int tbytes = TIFFTileSize(tif);
		t->data = xmalloc(tbytes);
		memset(t->data, 0, tbytes);
		int r = TIFFReadTile(tif, t->data, ii, jj, 0, 0);
		if (r != tbytes) fail("could not read tile");
	} else { // not tiled, read the whole image
		fail("image \"%s\" has no tiles!", filename);
		//int scanline_size = t->w * t->spp * t->bps/8;
		//for (int i = 0; i < t->h; i++)
		//	TIFFR
	}

	TIFFClose(tif);
}

// overwrite tile "idx" of the given file
static void insert_tile_into_file(char *filename, struct tiff_tile *t, int idx)
{
}

static void tiffu_whatever(char *filename)
{
	TIFF *tif = TIFFOpen(filename, "r");

	int w = tiff_imagewidth(tif);
	int h = tiff_imagelength(tif);
	int spp = tiff_samplesperpixel(tif);
	int Bps = tiff_bitspersample(tif)/8;
	int fmt = tiff_sampleformat(tif);

	if (TIFFIsTiled(tif)) {
		int tw = tiff_tilewidth(tif);
		int th = tiff_tilelength(tif);
		int ta = tiff_tilesacross(tif);
		int td = tiff_tilesdown(tif);
		int nt = TIFFNumberOfTiles(tif);
		bool broken = tiff_brokenpixels(tif);

		printf("TIFF w = %d\n", w);
		printf("TIFF h = %d\n", h);
		printf("TIFF spp = %d\n", spp);
		printf("TIFF Bps = %d\n", Bps);
		printf("TIFF fmt = %d\n", fmt);
		printf("TIFF broken pixels = %s\n", broken ? "yes" : "no");

		printf("TIFF tw = %d\n", tw);
		printf("TIFF th = %d\n", th);
		printf("TIFF ta = %d\n", ta);
		printf("TIFF td = %d\n", td);
		printf("TIFF ta * td = %d\n", ta * td);
		printf("TIFF nt = %d\n", nt);

		for (int j = 0; j < td; j++)
		for (int i = 0; i < ta; i++)
		{
			int ii = i*tw;
			int jj = j*th;
			int tidx = TIFFComputeTile(tif, ii, jj, 0, 0);
			fprintf(stderr, "%d %d : %d %d : %d (%d)\n",
					i, j, ii, jj, tidx, j*ta+i);
		}

		for (int t = 0; t < nt; t++)
		{
			int i = t % ta;
			int j = t / ta;
			int ii = i*tw;
			int jj = j*th;
			int tidx = TIFFComputeTile(tif, ii, jj, 0, 0);
			//fprintf(stderr, "t=%d i=%d j=%d tidx=%d\n",
			//		t, i, j, tidx);
			fprintf(stderr, "%d %d\n", t, tidx);
		}
	} else {
		printf("not tiled\n");
	}

	TIFFClose(tif);
}

static void tiffu_print_info(char *filename)
{
	TIFF *tif = TIFFOpen(filename, "r");
	if (!tif) fail("could not open TIFF file \"%s\"", filename);

	uint32_t w, h;
	uint16_t spp, bps, fmt;
	int r = 0;
	r += TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	r += TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	if (r != 2) fail("can not treat TIFF of unkwnown size");
	printf("TIFF %dx%d\n", w, h);

	r = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	if(!r) spp=1;
	if(r) printf("TIFF spp %d (r=%d)\n", spp, r);

	r = TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
	if(!r) bps=1;
	if(r) printf("TIFF bps %d (r=%d)\n", spp, r);

	r = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
	if(!r) fmt = SAMPLEFORMAT_UINT;
	if(r) printf("TIFF fmt %d (r=%d)\n", spp, r);

	if (TIFFIsTiled(tif)) {
		int tisize = TIFFTileSize(tif);
		uint32_t tilewidth, tilelength;
		r = 0;
		r += TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tilewidth);
		r += TIFFGetField(tif, TIFFTAG_TILELENGTH, &tilelength);
		if (r != 2) fail("TIFF tiles of unknown size");

		int ntiles = (w/tilewidth)*(h/tilelength);
		int ntiles2 = TIFFNumberOfTiles(tif);

		printf("TIFF ntiles %d (%d whole)\n", ntiles2, ntiles);
		printf("TIFF tilesize %d (%dx%d)\n", tisize,
				tilewidth, tilelength);

		int tiles_across = tiff_tilesacross(tif);
		int tiles_down = tiff_tilesdown(tif);
		int ntiles3 = tiles_across * tiles_down;

		printf("TIFF tile config = %dx%d (%d)\n",
				tiles_across, tiles_down, ntiles3);

	}

	TIFFClose(tif);
}

static int tiffu_ntiles(char *filename)
{
	int r = 1;

	TIFF *tif = TIFFOpen(filename, "r");
	if (!tif) fail("could not open TIFF file \"%s\"", filename);

	if (TIFFIsTiled(tif)) {
		uint32_t w, h, tw, tl;
		r += TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		r += TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		r += TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
		r += TIFFGetField(tif, TIFFTAG_TILELENGTH, &tl);
		if (r != 5) fail("TIFF tiles of unknown size");

		r = (w/tw)*(h/tl);
	}

	TIFFClose(tif);
	return r;
}

static int tiffu_ntiles2(char *filename)
{
	int r = 1;

	TIFF *tif = TIFFOpen(filename, "r");
	if (!tif) fail("could not open TIFF file \"%s\"", filename);

	//if (TIFFIsTiled(tif)) {
		r = TIFFNumberOfTiles(tif);
	//}

	TIFFClose(tif);
	return r;
}

// get a tile (low level: output image is filled to an array)
static float *tiffu_tget_ll(char *filename_in, int tile_idx,
		int *out_w, int *out_h, int *out_pd)
{

	TIFF *tif = TIFFOpen(filename_in, "r");
	if (!tif) fail("could not open TIFF file \"%s\"", filename_in);

	uint32_t w, h;
	uint16_t spp, bps, fmt;
	int r = 0;
	r += TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	r += TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	if (r != 2) fail("can not treat TIFF of unkwnown size");
	printf("TIFF %dx%d\n", w, h);

	r = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	if(!r) spp=1;
	if(r) printf("TIFF spp %d (r=%d)\n", spp, r);

	r = TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
	if(!r) bps=1;
	if(r) printf("TIFF bps %d (r=%d)\n", spp, r);

	r = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
	if(!r) fmt = SAMPLEFORMAT_UINT;
	if(r) printf("TIFF fmt %d (r=%d)\n", spp, r);

	if (TIFFIsTiled(tif)) {
		int tisize = TIFFTileSize(tif);
		uint32_t tilewidth, tilelength;
		r = 0;
		r += TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tilewidth);
		r += TIFFGetField(tif, TIFFTAG_TILELENGTH, &tilelength);
		if (r != 2) fail("TIFF tiles of unknown size");

		int ntiles = (w/tilewidth)*(h/tilelength);

		printf("TIFF ntiles %d\n", ntiles);
		printf("TIFF tilesize %d (%dx%d)\n", tisize,
				tilewidth, tilelength);

	}

	TIFFClose(tif);
	return NULL;
}

static void debug_tile_struct(struct tiff_tile *t)
{
	fprintf(stderr, "tile struct at %p\n", (void*)t);
	fprintf(stderr, "\tw = %d\n", t->w);
	fprintf(stderr, "\th = %d\n", t->h);
	fprintf(stderr, "\tspp = %d\n", t->spp);
	fprintf(stderr, "\tbps = %d\n", t->bps);
	fprintf(stderr, "\tfmt = %d\n", t->fmt);
}

// get a tile (high level: output image is saved to a file)
static void tiffu_tget_hl(char *filename_out, char *filename_in, int tile_idx)
{
	struct tiff_tile t[1];
	read_tile_from_file(t, filename_in, tile_idx);

	debug_tile_struct(t);

	write_tile_to_file(filename_out, t);
	free(t->data);
}

void tcrop(char *fname_out, char *fname_in, int x0, int xf, int y0, int yf)
{
	// open input file
	TIFF *tif = TIFFOpen(fname_in, "r");
	if (!tif) fail("can not open TIFF file \"%s\" for reading", fname_in);

	fprintf(stderr, "tif = %p\n", (void*)tif);

	// read input file info
	//struct tiff_file_info tinfo[1];
	struct tiff_file_info tinfo[1];
	get_tiff_file_info(tinfo, tif);


	// create output structure
	struct tiff_tile tout[1];
	tout->w = 1 + xf - x0;
	tout->h = 1 + yf - y0;
	tout->spp = tinfo->spp;
	tout->bps = tinfo->bps;
	tout->fmt = tinfo->fmt;
	int pixel_size = tinfo->spp * tinfo->bps/8;
	int output_size = tout->w * tout->h * pixel_size;
	//tout->data = xmalloc(tout->w * tout->h * pixel_size);
	tout->data = xmalloc(output_size);
	tout->broken = false;

	// convenience simplifcations (not necessary)
	if (!tinfo->tiled) fail("can only crop tiled files");
	if (tinfo->packed) fail("can not crop packed images");
	if (tinfo->broken) fail("can not crop images with broken pixels");

	// compute coordinates of required tiles
	int tix0 = x0 / tinfo->tw;
	int tixf = xf / tinfo->tw;
	int tiy0 = y0 / tinfo->th;
	int tiyf = yf / tinfo->th;
	assert(tix0 <= tixf); assert(tixf < tinfo->ta);
	assert(tiy0 <= tiyf); assert(tiyf < tinfo->td);

	// read input tiles
	int tisize = TIFFTileSize(tif);
	assert(tisize == tinfo->tw * tinfo->th * tinfo->spp * tinfo->bps/8);
	uint8_t *buf = xmalloc(tisize);

	// the outer 2 loops traverse the needed tiles of the big
	// the inner 3 loops traverse the needed pixels of each tile
	// i=tile indexes, ii=tile-wise pixel coords, iii=crop-wise pixel coords
	for (int j = tiy0; j <= tiyf; j++)
	for (int i = tix0; i <= tixf; i++)
	{
		TIFFReadTile(tif, buf, i*tinfo->tw, j*tinfo->th, 0, 0);
		int i_0 = i > tix0 ? 0 : x0 % tinfo->tw; // first pixel in tile
		int j_0 = j > tiy0 ? 0 : y0 % tinfo->th;
		int i_f = i < tixf ? tinfo->tw-1 : xf % tinfo->tw; // last pixel
		int j_f = j < tiyf ? tinfo->th-1 : yf % tinfo->th;
		assert(0 <= i_0); assert(i_0 < tinfo->tw);
		assert(0 <= j_0); assert(j_0 < tinfo->th);
		assert(0 <= i_f); assert(i_f < tinfo->tw);
		assert(0 <= j_f); assert(j_f < tinfo->th);
		for (int jj = j_0; jj <= j_f; jj++)
		for (int ii = i_0; ii <= i_f; ii++)
		for (int l = 0; l < pixel_size; l++)
		{
			int iii = ii + i*tinfo->tw - x0;
			int jjj = jj + j*tinfo->th - y0;
			int oidx = l + pixel_size*(iii + jjj*tout->w);
			int iidx = l + pixel_size*(ii + jj*tinfo->tw);
			assert(iii >= 0); assert(iii < tout->w);
			assert(jjj >= 0); assert(jjj < tout->h);
			assert(0 <= oidx); assert(oidx < output_size);
			assert(0 <= iidx); assert(iidx < tisize);
			tout->data[oidx] = buf[iidx];
		}
	}
	free(buf);

	// close input file
	TIFFClose(tif);

	// write output data
	write_tile_to_file(fname_out, tout);
	
	// cleanup
	free(tout->data);
}


//static void tiffu_tilesize(char *fname, int tile_idx, int *w, int *h, int *pd)
//{
//	TIFFSetWarningHandler(NULL);//suppress warnings
//
//	TIFF *tif = TIFFOpen(filename, "r");
//	if (!tif) fail("could not open TIFF file \"%s\"", filename);
//
//	uint32_t w, h;
//	uint16_t spp, bps, fmt;
//	int r = 0;
//	r += TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
//	r += TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
//	if (r != 2) fail("can not treat TIFF of unkwnown size");
//	printf("TIFF %dx%d\n", w, h);
//
//	r = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
//	if(!r) spp=1;
//	if(r) printf("TIFF spp %d (r=%d)\n", spp, r);
//
//	r = TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
//	if(!r) bps=1;
//	if(r) printf("TIFF bps %d (r=%d)\n", spp, r);
//
//	r = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
//	if(!r) fmt = SAMPLEFORMAT_UINT;
//	if(r) printf("TIFF fmt %d (r=%d)\n", spp, r);
//
//	if (TIFFIsTiled(tif)) {
//		int tisize = TIFFTileSize(tif);
//		uint32_t tilewidth, tilelength;
//		r = 0;
//		r += TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tilewidth);
//		r += TIFFGetField(tif, TIFFTAG_TILELENGTH, &tilelength);
//		if (r != 2) fail("TIFF tiles of unknown size");
//
//		int ntiles = (w/tilewidth)*(h/tilelength);
//
//		printf("TIFF ntiles %d\n", ntiles);
//		printf("TIFF tilesize %d (%dx%d)\n", tisize,
//				tilewidth, tilelength);
//
//	}
//
//	TIFFClose(tif);
//}

//static int tiffu_ntiles_cleant(char *filename)
//{
//	int r = 0;
//
//	TIFF *tif = TIFFOpen(filename, "r");
//
//	if (TIFFIsTiled(tif)) {
//		uint32_t tilewidth, tilelength;
//		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tilewidth);
//		TIFFGetField(tif, TIFFTAG_TILELENGTH, &tilelength);
//		r = (w/tilewidth)*(h/tilelength);
//	}
//
//	TIFFClose(tif);
//	return r;
//}

static int main_info(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s file.tiff\n", *v);
		return 1;
	}
	char *filename = v[1];

	tiffu_print_info(filename);

	return 0;
}

static void tiffu_print_info2(char *filename)
{
	struct tiff_file_info t[1];

}

static int main_info2(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s file.tiff\n", *v);
		return 1;
	}
	char *filename = v[1];

	tiffu_print_info2(filename);

	return 0;
}

static int main_whatever(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s file.tiff\n", *v);
		return 1;
	}
	char *filename = v[1];

	tiffu_whatever(filename);

	return 0;
}

static int main_ntiles(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s file.tiff\n", *v);
		return 1;
	}

	printf("%d (%d)\n", tiffu_ntiles(v[1]), tiffu_ntiles2(v[1]));
	//printf("%d\n", tiffu_ntiles2(v[1]));

	return 0;
}

static int main_tget(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s file.tiff idx til.tiff\n", *v);
		//                          0 1         2   3
		return 1;
	}
	char *filename_in = v[1];
	int tile_idx = atoi(v[2]);
	char *filename_out = v[3];

	tiffu_tget_hl(filename_out, filename_in, tile_idx);

	return 0;
}

static int main_crop(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s cx cy r in.tiff out.tiff\n", *v);
		//                          0 1  2  3 4       5
		return 1;
	}
	int cx = atoi(v[1]);
	int cy = atoi(v[2]);
	int rad = atoi(v[3]);
	char *filename_in = v[4];
	char *filename_out = v[5];

	int xmin = cx - rad; if (xmin < 0) xmin = 0;
	int ymin = cy - rad; if (ymin < 0) ymin = 0;
	int xmax = cx + rad; //if (xmax >= t->w) xmax = t->w - 1;
	int ymax = cy + rad; //if (ymax >= t->h) ymax = t->h - 1;

	tcrop(filename_out, filename_in, xmin, xmax, ymin, ymax);

	return 0;
}

static void my_putchar(FILE *f, int c)
{
	fputc(c, f);
}

static int main_imprintf(int argc, char *argv[])
{
	if (argc != 3) {
		fprintf(stderr, "usage:\n\t%s format image\n", *argv);
		//                          0 1      2
		return 1;
	}
	char *format = argv[1];
	char *filename = argv[2];

	struct tiff_file_info t;
	get_tiff_file_info_filename(&t, filename);

	FILE *f = stdout;
	while (*format) {
		int c = *format++;
		if (c == '%') {
			c = *format++; if (!c) break;
			int p = -1;
			switch(c) {
			case 'w': p = t.w; break;
			case 'h': p = t.h; break;
			case 'n': p = t.w * t.h; break;
			case 'd': p = t.spp; break;
			case 'b': p = t.bps; break;
			case 'B': p = t.bps/8; break;
			case 'p': p = t.spp * t.bps; break;
			case 'P': p = t.spp * t.bps / 8; break;
			case 'f': p = t.fmt == SAMPLEFORMAT_IEEEFP; break;
			case 'u': p = t.fmt == SAMPLEFORMAT_UINT; break;
			case 'i': p = t.fmt == SAMPLEFORMAT_INT; break;
			case 'T': p = t.tiled; break;
			case 'Z': p = t.broken; break;
			case 'W': p = t.tw; break;
			case 'H': p = t.th; break;
			case 'A': p = t.ta; break;
			case 'D': p = t.td; break;
			default: break;
			}
			fprintf(f, "%d", p);
		} else if (c == '\\') {
			c = *format++; if (!c) break;
			if (c == '%') my_putchar(f, '%');
			if (c == 'n') my_putchar(f, '\n');
			if (c == 't') my_putchar(f, '\t');
			if (c == '"') my_putchar(f, c);
			if (c == '\'') my_putchar(f, c);
			if (c == '\\') my_putchar(f, c);
		} else
			my_putchar(f, c);
	}

	return 0;
}

int main(int c, char *v[])
{
	TIFFSetWarningHandler(NULL);//suppress warnings

	if (c < 2) {
	err:	fprintf(stderr, "usage:\n\t%s {info|ntiles|tget|tput}\n", *v);
		return 1;
	}

	if (0 == strcmp(v[1], "info"))     return main_info(c-1, v+1);
	if (0 == strcmp(v[1], "imprintf")) return main_imprintf(c-1, v+1);
	if (0 == strcmp(v[1], "ntiles"))   return main_ntiles(c-1, v+1);
	if (0 == strcmp(v[1], "tget"))     return main_tget(c-1, v+1);
	if (0 == strcmp(v[1], "whatever")) return main_whatever(c-1, v+1);
	if (0 == strcmp(v[1], "crop"))     return main_crop(c-1, v+1);
	//if (0 == strcmp(v[1], "tput"))   return main_tput(c-1, v+1);

	goto err;
}
