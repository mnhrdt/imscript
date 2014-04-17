// tiff utils
//
// info   f.tiff
// ntiles f.tiff
// tget   f.tiff n til.png
// tput   f.tiff n til.png


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>

#include "fail.c"

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

		printf("TIFF ntiles %d\n", ntiles);
		printf("TIFF tilesize %d (%dx%d)\n", tisize,
				tilewidth, tilelength);

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

#include "iio.h"

// get a tile (high level: output image is saved to a file)
static void tiffu_tget_hl(char *filename_out, char *filename_in, int tile_idx)
{
	int w, h, pd;
	float *tile_data = tiffu_tget_ll(filename_in, tile_idx, &w, &h, &pd);
	iio_save_image_float_vec(filename_out, tile_data, w, h, pd);
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

static int main_ntiles(int c, char *v[])
{
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s file.tiff\n", *v);
		return 1;
	}

	printf("%d\n", tiffu_ntiles(v[1]));

	return 0;
}

static int main_tget(int c, char *v[])
{
	if (c != 4) {
		fprintf(stderr, "usage:\n\t%s file.tiff idx til.png\n", *v);
		//                          0 1         2   3
		return 1;
	}
	char *filename_in = v[1];
	int tile_idx = atoi(v[2]);
	char *filename_out = v[3];

	int ntiles = tiffu_ntiles(filename_in);
	if (tile_idx < 0 || tile_idx >= ntiles)
		tile_idx = 0;

	tiffu_tget_hl(filename_in, tile_idx, filename_out);

	printf("%d\n", tiffu_ntiles(v[1]));

	return 0;
}

int main(int c, char *v[])
{
	if (c < 2) {
	err:	fprintf(stderr, "usage:\n\t%s {info|ntiles|tget|tput}\n", *v);
		return 1;
	}
	TIFFSetWarningHandler(NULL);//suppress warnings
	if (0 == strcmp(v[1], "info"))   return main_info(c-1, v+1);
	if (0 == strcmp(v[1], "ntiles")) return main_ntiles(c-1, v+1);
	if (0 == strcmp(v[1], "tget"))   return main_tget(c-1, v+1);
	//if (0 == strcmp(v[1], "tput"))   return main_tput(c-1, v+1);
	goto err;
}
