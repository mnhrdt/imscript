#!/bin/sh



# parameters for the matching
M_MATCH=290
M_RANGE=200

# parameters for ransac
RANSAC_N=10000
RANSAC_E=0.5
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
GEOMETRY=`imprintf '%w %h' $FIRSTIMAGE`

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
		graysift.sh lowe $PGMFILE > $SIFTFILE
	fi
}

# STEP 1: compute the (possibly cached) sift keypoints of each image
for i in `seq $FIRST $LAST`; do
	INPNG=`printf $INPAT $i`
	OUTPNG=`printf $OUTPAT $i`
	INSIFT=`printf $SIFTPAT $i`.sift
	compute_sift_from_png $INPNG $INSIFT
done

# STEP 2: match the keypoints of each image to those of the first one
for i in `seq $FIRST $LAST`; do
	SIFT0=`printf $SIFTPAT $FIRST`.sift
	SIFTi=`printf $SIFTPAT $i`.sift
	PAIRS0i=`printf $REGPAT $i`.pairs
	VPAIRS0i=view_pairs_`printf $REGPAT $i`.png
	echo "matching image \"$i\" to image \"$FIRST\""
	if [ "$i" = "$FIRST" ] ; then
		#echo "\tthey are the same, do nothing here"
		true
	else
		if [ ! -s $PAIRS0i ]; then
		srmatch $M_MATCH $SIFT0 $SIFTi 200 2 $PAIRS0i /dev/null /dev/null
		pview pairs 1 0 0 0 1 0 0 0 1 $GEOMETRY < $PAIRS0i > $VPAIRS0i
		#siftu pairr $M_MATCH $SIFT0 $SIFTi $M_RANGE $M_RANGE $PAIRS0i
		fi
	fi
done

# STEP 3: find homographies among the pairs of keypoints above
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	SIFT0=`printf $SIFTPAT $FIRST`.sift
	SIFTi=`printf $SIFTPAT $i`.sift
	PAIRS0i=`printf $REGPAT $i`.pairs
	HOM0i=`printf $REGPAT $i`.hom
	OMASK0i=`printf $REGPAT $i`.omask
	VOPAIRS0i=view_opairs_`printf $REGPAT $i`.png
	RR0i=`printf $REGPAT $i`.png
	echo "matching image \"$i\" to image \"$FIRST\""
	if [ "$i" = "$FIRST" ] ; then
		#echo "\tthey are the same, do nothing here"
		echo "1 0 0 0 1 0 0 0 1" > $HOM0i
		cp $PNGi $RR0i
	else
		if [ ! -s $HOM0i ]; then
		ransac hom $RANSAC_N $RANSAC_E $RANSAC_M < $PAIRS0i
		cp /tmp/omask.txt $OMASK0i
		cp /tmp/ramo.txt $HOM0i
		pview pairs 1 0 0 0 1 0 0 0 1 $GEOMETRY $OMASK0i < $PAIRS0i > $VOPAIRS0i
		synflow homi "`cat $HOM0i`" $PNGi $RR0i /dev/null
		fi
	fi
done

