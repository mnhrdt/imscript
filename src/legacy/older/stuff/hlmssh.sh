#!/bin/bash

# parameters for low-resolution initial match
SMTRES=290
TOP=600
RAD=6

# parameters for full-resolution refinement
NTRIALS=10000
MAXERR=0.5
MINLIERS=1000

set -e

if [ "$#" != "3" ]; then
	echo -e "usage:\n\t$0 "\
	"a.png b.png {hom|aff} > model.txt"
	# 1    2      3
	false
fi

IMGA=$1
IMGB=$2
MODEL=$3

T=`mktemp -d /tmp/hlmssh.XXXXXX` || exit 1
echo "TPD=$T" 1>&2
echo "$@" > $T/args

fnorm $IMGA | qauto - $T/a.pgm
fnorm $IMGB | qauto - $T/b.pgm
graysift.sh zero $T/a.pgm > $T/a.sift
graysift.sh zero $T/b.pgm > $T/b.sift

srmatch $SMTRES $T/a.sift $T/b.sift $TOP $RAD $T/p.txt $T/m.txt $T/h.txt
ransac hom $NTRIALS $MAXERR $MINLIERS $T/rh.txt $T/rm.txt $T/ri.txt < $T/p.txt

cat $T/rh.txt
