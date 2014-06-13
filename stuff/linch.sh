#!/bin/bash

# linear contrast enhancement

set -e

if [ "$#" != "4" ]; then
	echo -e "usage:\n\t$0 "\
	"lambda sigma prec in > out"
	# 1     2     3    4
	false
fi

LAMBDA=$1
SIGMA=$2
PREC=$3
IMGIN=$4
export FIWARN=0

PREC2=`plambda -c $PREC 'dup *'`
LAMBDA2=`plambda -c $LAMBDA 'dup *'`

TPD=`mktemp -d /tmp/linch.XXXXXX` || exit 1
#echo "TPD=$TPD" 1>&2
#echo "$@" > $TPD/args

fft < $IMGIN > $TPD/imgin.fft
plambda $IMGIN ":i :j hypot 0.1 <" | blur l $SIGMA | fft | plambda - $TPD/imgin.fft "x[0] dup dup * $PREC2 + / :I :I * :J :J * + dup $LAMBDA2 + / y * *" 2>/dev/null | ifft

