// SRTM4 files are of size 6000x6000 pixels and cover an area of 5x5 degrees

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tiffio.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef TIFFU_C_INCLUDED
#define TIFFU_OMIT_MAIN
#include "tiffu.c"
#endif//TIFFU_C_INCLUDED

#define SRTM4_TILE_MEGABYTES 2


#define NO_DATA 0

#define SRTM4_TIF "%s/srtm_%02d_%02d.tif"
#define SRTM4_TIFO "%s/srtm_%02d_%02d.tif.%%d"

//#define SRTM4_URL_TIF "http://138.231.80.250:443/srtm/tiff/srtm_%02d_%02d.zip"
 #define SRTM4_URL_TIF "ftp://xftp.jrc.it/pub/srtmV4/tiff/srtm_%02d_%02d.zip"
// #define SRTM4_URL_TIF "--http-user=data_public --http-password=GDdci http://data.cgiar-csi.org/srtm/tiles/GeoTIFF/srtm_%02d_%02d.zip"



// download the contents of an url into a file
static int download(const char *to_file, const char *from_url)
{
	int nbuf = 2*FILENAME_MAX;
	char buf[nbuf];
	snprintf(buf, nbuf, "wget %s -O %s", from_url, to_file);
	int r = system(buf);
	return r;
}

// return the name of the cache directory
static char *cachedir(void)
{
	// XXX TODO FIXME WRONG WARNING : THIS code is not reentrant
	static char output_dirname[FILENAME_MAX];
	char *env_cache = getenv("SRTM4_CACHE");
	if (env_cache) {
		snprintf(output_dirname, FILENAME_MAX, "%s", env_cache);
	} else {
		char *homedir = getenv("HOME");
		if (!homedir)
			homedir = "/tmp";
		snprintf(output_dirname, FILENAME_MAX, "%s/.srtm4", homedir);
	}
	int r = mkdir(output_dirname, 0777);
	if (r != 0 && errno != EEXIST) {
		fprintf(stderr, "SRTM4: cannot create cachedir \"%s\"\n",
				output_dirname);
		exit(2);
	}
	return output_dirname;
}


static void get_tile_index_and_position(
		int *tlon, int *tlat,
		float *xlon, float *xlat,
		double lon, double lat)
{
	if (lat > 60) lat = 60;
	if (lat < -60) lat = -60;

	// tiles longitude indexes go from 1 to 72, covering the range from
	// -180 to +180
	*tlon = fmod(1+floor((lon+180)/5), 72);
	if (*tlon == 0) *tlon = 72;
	lon = fmod(lon + 180, 360);
	*xlon = 1200*fmod(lon, 5);

	// tiles latitude indexes go from 1 to 24, covering the range from 60
	// to -60
	*tlat = 1+floor((60-lat)/5);
	if (*tlat == 25) *tlat = 24;
	*xlat = 1200*fmod(60-lat, 5);

//	fprintf(stderr, "tlon = %d\n", *tlon);
//	fprintf(stderr, "tlat = %d\n", *tlat);
//	fprintf(stderr, "xlon = %g\n", *xlon);
//	fprintf(stderr, "xlat = %g\n", *xlat);

	assert(*tlon > 0); assert(*tlon <= 72);
	assert(*tlat > 0); assert(*tlat <= 24);
	assert(*xlon >= 0); assert(*xlon < 6000);
	assert(*xlat >= 0); assert(*xlat < 6000);
}

static bool file_exists(const char *fname)
{
	FILE *f = fopen(fname, "r");
	if (f)
	{
		fclose(f);
		return true;
	}
	return false;
}

static void download_tile_file(int tlon, int tlat)
{
	int n = FILENAME_MAX;
	char url[n], zipname[n];

	snprintf(url, n, SRTM4_URL_TIF, tlon, tlat);
	snprintf(zipname, n, "%s/tmp.zip", cachedir());
	int rd = download(zipname, url);
	if (0 == rd) {
		char cmd[n];
		snprintf(cmd, n, "unzip -qq -o %s -d %s", zipname, cachedir());
		rd = system(cmd);
		if (rd) {
			fprintf(stderr,"failed unzipping file %s\n", zipname);
		}
		// TODO: do something if unzip fails
	}
}

static char *get_tile_filename(int tlon, int tlat, bool check_octaves)
{
	// XXX TODO FIXME WRONG WARNING : THIS code is not reentrant
	static char fname[FILENAME_MAX];
	char *fmt = check_octaves ? SRTM4_TIFO : SRTM4_TIF;
	snprintf(fname, FILENAME_MAX, fmt, cachedir(), tlon, tlat);
	return fname;
}

static struct tiff_octaves *global_table_of_tiles[360][180] = {{0}};
static int global_table_of_count[360][180] = {{0}};

static struct tiff_octaves *produce_tile(int tlon, int tlat)
{
	struct tiff_octaves *t = global_table_of_tiles[tlon][tlat];
	if (!t) {
		char *fname = get_tile_filename(tlon, tlat, true);
		char fname0[FILENAME_MAX];
		snprintf(fname0, FILENAME_MAX, fname, 0);
		fprintf(stderr, "trying file \"%s\"\n", fname0);
		if (!file_exists(fname0)) {
			fname = get_tile_filename(tlon, tlat, false);
			fprintf(stderr, "trying now file \"%s\"\n", fname);
			if (!file_exists(fname)
				&& global_table_of_count[tlon][tlat] < 2)
			{
				global_table_of_count[tlon][tlat] += 1;
				download_tile_file(tlon, tlat);
			}
		}
		if (!file_exists(fname) && !file_exists(fname0))
		{
			//fprintf(stderr, "WARNING: srtm4 tile \"%d %d\" "
			//		"not available\n", tlon, tlat);
			return NULL;
		}
		t = malloc(sizeof*t);
		tiff_octaves_init(t, fname, 0);
		if ((t->i->w != 6000) || (t->i->h != 6000))
			exit(fprintf(stderr, "produce_tile: srtm4 file not "
					"6000x6000\n"));
		global_table_of_tiles[tlon][tlat] = t;
	}
	return t;
}

static double getpixelo_double(struct tiff_octaves *t, double x, double y, int o)
{
	if (o > t->noctaves)
		o = t->noctaves - 1;
	float ofac = 1 << o;
	int i = x / ofac;
	int j = y / ofac;
	assert(t->i->fmt == SAMPLEFORMAT_INT);
	assert(t->i->bps == 16);
	assert(t->i->spp == 1);
	int16_t *pixel = tiff_octaves_getpixel(t, o, i, j);
	return (double)*pixel;
}

double srtm4o(double lon, double lat, int octave)
{
	if (lat >= 60 || lat <= -60)
		return NO_DATA;
	int tlon, tlat;
	float xlon, xlat;
	get_tile_index_and_position(&tlon, &tlat, &xlon, &xlat, lon, lat);
	if (octave < 0)
		fprintf(stderr, "lonlat: %g %g  => xlonxlat: %g %g\n", lon, lat, xlon, xlat);

	struct tiff_octaves *t = produce_tile(tlon, tlat);
	if (t == NULL)
		return NO_DATA;
	else {
		double r = getpixelo_double(t, xlon, xlat, fabs(octave));
		return r > 0 ? r : 0;
	}
}

void srtm4_free_tiles(void)
{
	for (int j = 0; j < 360; j++)
	for (int i = 0; i < 180; i++)
	if (global_table_of_tiles[j][i])
	{
		tiff_octaves_free(global_table_of_tiles[j][i]);
			free(global_table_of_tiles[j][i]);
	}
}

#ifdef MAIN_SRTM4
int main(int c, char *v[])
{
	if (c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s longitude latitude\n", *v);
		return 1;
	}
    if (c == 3) {
	    double lon = atof(v[1]);
	    double lat = atof(v[2]);
	    double r = srtm4_wrt_ellipsoid(lon, lat);
	    printf("%g\n", r);
	    return 0;
    }
    else {
        double lon, lat, r;
        while(2 == scanf("%lf %lf\n", &lon, &lat)) {
            r = srtm4_nn_wrt_ellipsoid(lon, lat);
            printf("%g\n", r);
        }
    }
}
#endif//MAIN_SRTM4


#ifndef MAIN_SRTM4
#ifdef MAIN_SRTM4_WHICH_TILE
static void print_tile_filename(double lon, double lat)
{
	int tlon, tlat;
	float xlon, xlat;
	get_tile_index_and_position(&tlon, &tlat, &xlon, &xlat, lon, lat);
	printf("srtm_%02d_%02d\n", tlon, tlat);
	return;
}


int main(int c, char *v[])
{
	if (c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s longitude latitude\n", *v);
		return 1;
	}
    if (c == 3) {
	    double lon = atof(v[1]);
	    double lat = atof(v[2]);
        print_tile_filename(lon, lat);
	    return 0;
    }
    else {
        double lon, lat;
        while(2 == scanf("%lf %lf\n", &lon, &lat)) {
            print_tile_filename(lon, lat);
        }
    }
}
#endif//MAIN_SRTM4_WHICH_TILE
#endif
