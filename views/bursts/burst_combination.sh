#!/bin/sh


# parameters for computing the wheights
W_WINSIZE=99
W_TOPQ=30

# stencil for computing the local contrast measure
NGRAD="x(1,0) x - x(0,1) x - hypot split + + 3 /"


set -e

if [ "$#" != "4" ]; then
	echo "usage:\n\t$0 inpat first last outimg"
	exit 1
fi

INPAT=$1
FIRST=$2
LAST=$3
OUTIMG=$4

echo "INPAT = $INPAT"
echo "FIRST = $FIRST"
echo "LAST = $LAST"
echo "OUTIMG = $OUTIMG"

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



# required infrastructure from IMSCRIPT
testex plambda
testex fft
testex ifft
testex qeasy
testex vecov
testex blur



# common filename patterns
#SIFTPAT=`basename $INPAT .png`
#REGPAT=`basename $OUTPAT .png`
#FIRSTIMAGE=`printf $INPAT $FIRST`
#GEOMETRY=`imprintf '%w %h' $FIRSTIMAGE`

# print the file names for debugging
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	wPNGi=wD_`printf $INPAT $i`
	echo "$i: $PNGi $wPNGi"
done

# compute the weights of each image
ALL_I=""
ALL_W=""
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	wPNGi=wD_`printf $INPAT $i`
	ALL_I="$ALL_I $PNGi"
	ALL_W="$ALL_W $wPNGi"
	if [ ! -s $wPNGi ]; then
		echo "computing weights for image \"$PNGi\""
		plambda $PNGi "$NGRAD" | blur square $W_WINSIZE | qeasy 0 $W_TOPQ - $wPNGi
	fi
done

# pixelwise sum of all the weights
vecov sum $ALL_W > sum_wD

# pixelwise linear combinations of weight plus image
ALL_XW=""
for i in `seq $FIRST $LAST`; do
	PNGi=`printf $INPAT $i`
	wPNGi=wD_`printf $INPAT $i`
	xwPNGi=xwD_`printf $INPAT $i`
	ALL_XW="$ALL_XW $xwPNGi"
	if [ ! -s $xwPNGi ]; then
		plambda $PNGi $wPNGi "x y *" > $xwPNGi
	fi
done

# unnormalized sum of all the weighted images
vecov sum $ALL_XW > sum_xwD

# normalized sum of all the weighted images
plambda sum_xwD sum_wD "x y /" | qeasy 0 255 - $OUTIMG
