#!/bin/bash

# registration by sift + ransac homography

set -x
set -e

# sift matching parameters
SIFTDIST=290
PSTART=1000
IERR=7

# ransac parameters
NTRIALS=10000
RERR=1.5
MINMATCH=30


if [ "$#" != "3" ]; then
	echo -e "usage:\n\t$0 "\
	"a.png b.png reg_b.png"
	#1     2     3
	false
fi

TPD=`mktemp -d /tmp/homregister.XXXXXX` || exit 1
echo "TPD=$TPD" 1>&2
echo "$@" > $TPD/args

FILE_A=$1
FILE_B=$2
FILE_OUTPUT=$3


cp $FILE_A $TPD/a.png
cp $FILE_B $TPD/b.png

unisift.sh lowe $FILE_A > $TPD/a.sift
unisift.sh lowe $FILE_B > $TPD/b.sift
srmatch    $SIFTDIST $TPD/{a,b}.sift $PSTART $IERR $TPD/o{p,m,h}.txt
#siftu pairr $SIFTDIST $TPD/a.sift $TPD/b.sift $IERR $IERR $TPD/ab.pairs
ransac hom $NTRIALS $RERR $MINMATCH $TPD/r{h,m,i}.txt < $TPD/op.txt
synflow homi "`cat $TPD/rh.txt`" $TPD/{,r}b.png /dev/null

cp $TPD/rb.png $FILE_OUTPUT
