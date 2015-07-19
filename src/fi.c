#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fi.h"
#include "iio.h"
#define TIFFU_OMIT_MAIN
#include "tiff_octaves.c"

// the following struct is an implementation detail,
// it is only used on this file
//
struct FI {
	// part of the public interface, "struct fancy_image"
	int w;
	int h;
	int pd;
	int no;

	// implementation details
	bool tiffo;
	struct tiff_octaves t[1];
	float *x;
};

// Compiler trick to check that "FI" can fit inside a "fancy_image"
// (if this line fails, increase the padding at the struct fancy_image on fi.h)
typedef char check_FI_size[sizeof(struct FI)<=sizeof(struct fancy_image)?1:-1];


// check whether a filename corresponds to a small image or a tiled tiff
static bool filename_corresponds_to_tiffo(char *filename)
{
	if (0 == strcmp(filename, "-"))
		return false;
	struct tiff_info ti[1];
	bool r = get_tiff_info_filename_e(ti, filename);
	if (!r)
		return false;
	return ti->tiled;
}

// API: open a fancy image
struct fancy_image fancy_image_open(char *filename, double megabytes)
{
	struct fancy_image r[1];
	struct FI *f = (void*)r;
	if (filename_corresponds_to_tiffo(filename)) {
		f->tiffo = true;
		tiff_octaves_init(f->t, filename, megabytes);
		f->w = f->t->i->w;
		f->h = f->t->i->h;
		f->pd = f->t->i->spp;
		f->no = f->t->noctaves;
	} else {
		f->tiffo = false;
		f->x = iio_read_image_float_vec(filename, &f->w, &f->h, &f->pd);
	}
	return *r;
}

// API: close a fancy image
void fancy_image_close(struct fancy_image *fi)
{
	struct FI *f = (void*)fi;

	if (f->tiffo)
		tiff_octaves_free(f->t);
	else
		free(f->x);
}

// internal conversion function function
static float convert_sample_to_float(struct tiff_info *ti, void *p)
{
	switch(ti->fmt) {
	case SAMPLEFORMAT_UINT:
		if (8 == ti->bps)        return *(uint8_t *)p;
		else if (16 == ti->bps)  return *(uint16_t *)p;
		else if (32 == ti->bps)  return *(uint16_t *)p;
		break;
	case SAMPLEFORMAT_INT:
		if (8 == ti->bps)        return *(int8_t *)p;
		else if (16 == ti->bps)  return *(int16_t *)p;
		else if (32 == ti->bps)  return *(int16_t *)p;
		break;
	case SAMPLEFORMAT_IEEEFP:
		if (32 == ti->bps)       return *(float *)p;
		else if (64 == ti->bps)  return *(double *)p;
		break;
	}
	return NAN;
}

// API: get a sample from an image
float fancy_image_getsample_float(struct fancy_image *fi, int i, int j, int l)
{
	struct FI *f = (void*)fi;

	if (i < 0 || j < 0 || i >= f->w || j >= f->h)
		return NAN;
	if (l < 0) l = 0;
	if (l < f->pd) l = f->pd - 1;

	if (f->tiffo) {
		uint8_t *p_pixel = tiff_octaves_getpixel(f->t, 0, i, j);
		uint8_t *p_sample = p_pixel + (l * f->t->i->bps) / 8;
		return convert_sample_to_float(f->t->i, p_sample);
	} else {
		int idx = (j * f->w + i) * f->pd + l;
		return f->x[idx];
	}
}

#ifdef MAIN_FI
#include <stdio.h>
int main(int c, char *v[])
{
	// process input arguments
	if (c != 2) {
		fprintf(stderr, "usage:\n\t%s image\n", *v);
		return 1;
	}
	char *filename = v[1];

	// do stuff
	double megabytes = 100;
	struct fancy_image f = fancy_image_open(filename, megabytes);
	printf("image \"%s\" : %dx%d, %d\n", filename, f.w, f.h, f.pd);
	for (int x = 1; x < 100000; x *= 10)
	{
		float s = fancy_image_getsample_float(&f, x, x, 0);
		printf("\tsample(%d,%d,0) = %g\n", x, x, s);
	}
	fancy_image_close(&f);

	// exit
	return 0;
}
#endif
