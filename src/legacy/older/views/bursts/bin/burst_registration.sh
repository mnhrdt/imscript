#!/bin/sh



# parameters for the matching
M_MATCH=290
M_RANGE=200

SRM_FIRST=600
SRM_E=4

# parameters for ransac
RANSAC_N=1000
RANSAC_E=0.9
RANSAC_M=30


set -e

if [ "$#" != "4" ]; then
	echo "usage:\n\t$1 inpat outpat first last"
	exit 1
fi

INPAT=$1
OUTPAT=$2
FIRST=$3
LAST=$4

echo "INPAT = $INPAT"
echo "OUTPAT = $OUTPAT"
echo "FIRST = $FIRST"
echo "LAST = $LAST"

# utility function "testex"
testex() {
	if which $1 > /dev/null ; then
		echo > /dev/null
	else
		echo "ERROR: required executable \"$1\" not available" >&2
		exit 1
	fi
}

# required infrastructure from UNIX
testex seq
testex basename
testex printf

# required infrastructure for SIFT
testex graysift.sh

# required infrastructure from IMSCRIPT
testex siftu
testex ransac
testex synflow
testex imprintf









# common filename patterns
SIFTPAT=`basename $INPAT .png`
REGPAT=`basename $OUTPAT .png`
FIRSTIMAGE=`printf $INPAT $FIRST`
GEOFILE=.geometry

if [ ! -s $GEOFILE ]; then
	imprintf "%w %h\n" $FIRSTIMAGE > $GEOFILE
fi
GEOMETRY=`cat $GEOFILE`

# print the file names for debugging
for i in `seq $FIRST $LAST`; do
	INI=`printf $INPAT $i`
	OUTI=`printf $OUTPAT $i`
	SIFTI=`printf $SIFTPAT $i`.sift
	echo "$i: $INI $SIFTI $OUTI"
done

# ALGORITM "burst registration"
# 1. compute sift keypoints of all images (may use caching)
# 2. match the keypoints of each image to those of the first one
# 3. run ransac to find an homography among each collection of pairs
# 3. resample each image to the first one using the computed homographies

compute_sift_from_png() {
	PNGFILE=$1
	SIFTFILE=$2

	if [ -s $SIFTFILE ]; then
		echo "using cached \"$SIFTFILE\" as sift of image \"$PNGFILE\""
	else
		PGMFILE=`basename $PNGFILE .png`.pgm
		plambda $PNGFILE "x split + + 3 /" | qeasy 0 255 - $PGMFILE
		graysift.sh zero $PGMFILE > $SIFTFILE
	fi
}

# STEP 1: compute the (possibly cached) sift keypoints of each image
for i in `seq $FIRST $LAST`; do
	INPNG=`printf $INPAT $i`
	OUTPNG=`printf $OUTPAT $i`
	INSIFT=`printf $SIFTPAT $i`.sift
	compute_sift_from_png $INPNG $INSIFT #&
done
#wait

# STEP 2: match the keypoints of each image to those of the first one
for i in `seq $FIRST $LAST`; do
	SIFT0=`printf $SIFTPAT $FIRST`.sift
	SIFTi=`printf $SIFTPAT $i`.sift
	PAIRSi=`printf $REGPAT $i`.pairs
	VPAIRSi=view_pairs_`printf $REGPAT $i`.png
	echo "roughly matching image \"$i\" to image \"$FIRST\""
	if [ "$i" = "$FIRST" ] ; then
		#echo "\tthey are the same, do nothing here"
		true
	else
		if [ ! -s $PAIRSi ]; then
		srmatch $M_MATCH $SIFT0 $SIFTi $SRM_FIRST $SRM_E $PAIRSi /dev/null /dev/null &
		#pview pairs 1 0 0 0 1 0 0 0 1 $GEOMETRY < $PAIRSi > $VPAIRSi
		#siftu pairr $M_MATCH $SIFT0 $SIFTi $M_RANGE $M_RANGE $PAIRSi
		fi
	fi
done
wait

# STEP 3: find homographies among the pairs of keypoints above
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	SIFT0=`printf $SIFTPAT $FIRST`.sift
	SIFTi=`printf $SIFTPAT $i`.sift
	PAIRSi=`printf $REGPAT $i`.pairs
	HOMi=`printf $REGPAT $i`.hom
	MASKi=`printf $REGPAT $i`.omask
	VOPAIRSi=view_opairs_`printf $REGPAT $i`.png
	RRi=`printf $REGPAT $i`.png
	echo "matching image \"$i\" to image \"$FIRST\""
	if [ "$i" = "$FIRST" ] ; then
		#echo "\tthey are the same, do nothing here"
		echo "1 0 0 0 1 0 0 0 1" > $HOMi
		cp $PNGi $RRi
	else
		if [ ! -s $HOMi ]; then
		ransac hom $RANSAC_N $RANSAC_E $RANSAC_M $HOMi $MASKi < $PAIRSi &
		pview pairs 1 0 0 0 1 0 0 0 1 $GEOMETRY $MASKi < $PAIRSi > $VOPAIRSi &
		#echo "\tRESAMPLING"
		#synflow homi "`cat $HOMi`" $PNGi $RRi /dev/null #&
		fi
	fi
done
wait

# STEP 4: resample the images by the obtained homographies
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	HOMi=`printf $REGPAT $i`.hom
	RRi=`printf $REGPAT $i`.png
	echo "moving image \"$i\" to image \"$FIRST\""
	if [ "$i" = "$FIRST" ] ; then
		#echo "\tthey are the same, do nothing here"
		cp $PNGi $RRi
	else
		if [ ! -s $RRi ]; then
			synflow homi "`cat $HOMi`" $PNGi $RRi /dev/null &
		fi
	fi
done
wait
