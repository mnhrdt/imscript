#!/bin/bash

set -x
set -e

# tile side
OPT_N=256

# pyramid factors
OPT_L="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192"

# compression
OPT_K=LZW
OPT_P=2


# XXX TODO FIXME : this code uses gdaltranslate and gdaladdo
# we should rewrite a standalone version that does exactly the same thing

# extract input arguments
if [ $# != "2" ] ; then
	printf "usage:\n\t$0 in.tif out.tif\n" >&2
	#                    1      2
	exit 1
fi
FNAME_IN=$1
FNAME_OUT=$2

TPD=`mktemp -d` || exit 2
FNAME_TMP=$TPD/a.tif

TILE_OPTIONS="-co TILED=YES -co BLOCKXSIZE=$OPT_N -co BLOCKYSIZE=$OPT_N"

# 1. convert original image to tiles
gdal_translate $FNAME_IN  $FNAME_TMP  $TILE_OPTIONS

# 2. build pyramid
gdaladdo -r average $FNAME_TMP $OPT_L

# 3. retile the pyramid
PYR_OPTIONS="-co COPY_SRC_OVERVIEWS=YES --config GDAL_TIFF_OVR_BLOCKSIZE $OPT_N"
K_OPTIONS="-co COMPRESS=$OPT_K -co PREDICTOR=$OPT_P"
ALL_OPTIONS="$TILE_OPTIONS $PYR_OPTIONS $K_OPTIONS"
gdal_translate $FNAME_TMP $FNAME_OUT $ALL_OPTIONS


# cleanup and exit
rm $FNAME_TMP
rmdir $TPD
