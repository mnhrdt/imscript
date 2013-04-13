#!/bin/bash

set -e

if [ "$#" != "4" ]; then
	echo -e "usage:\n\t$0 "\
	"ktype[g|l|u|c] kvar prec blurry > sharp"
	# 1             2    3    4
	false
fi

KTYPE=$1
KVAR=$2
PREC=$3
BLURRY=$4
export FIWARN=0

TPD=`mktemp -d /tmp/dwiener.XXXXXX` || exit 1
#echo "TPD=$TPD" 1>&2
#echo "$@" > $TPD/args

fft < $BLURRY > $TPD/blurry.fft
plambda $BLURRY ":i :j hypot 0.1 <" | blur $KTYPE $KVAR | fft | plambda - $TPD/blurry.fft "x[0] dup dup * $PREC + / y *" 2>/dev/null | ifft

