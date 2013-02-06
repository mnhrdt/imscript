#!/bin/bash

set -e

if [ "$#" != "8" ]; then
	echo -e "usage:\n\t$0 "\
	"w h d graintype[g|l|u|c] grainsize dist[g|l|u|c] strength seed > out"
	#1 2 3 4                  5         6             7        8
	false
fi

WIDTH=$1
HEIGHT=$2
DIM=$3
GRAINTYPE=$4
GRAINSIZE=$5
DIST=$6
STRENGTH=$7
export SRAND=$8  # environment variable used by plambda
export FIWARN=0


if [ "$GRAINTYPE" != "g" ]; then
	echo -e "ERROR($0): only 'g' graintype is supported for now" >&2
	false
fi

if [ "$DIM" == "1" ]; then
	case "$DIST" in
	"g")	RANDSOURCE="randn"
		PLAMBDA_SECOND="" ;;
	"u")	RANDSOURCE="randn"
		PLAMBDA_SECOND="2 sqrt / erf 3 sqrt *" ;;
	"c")	RANDSOURCE="randc"
		PLAMBDA_SECOND="0.01 *" ;;
	"l")	RANDSOURCE="randn randn randn randn 4 njoin"
		PLAMBDA_SECOND="split * >1 * <1 - 2 sqrt /" ;;
	*)	echo -e "ERROR($0): unrecognized dist '$DIST'" >&2 false ;;
	esac
	PLAMBDA_FIRST="$RANDSOURCE pi sqrt * 2 * $GRAINSIZE *"
	PLAMBDA_SECOND="$PLAMBDA_SECOND $STRENGTH *"
elif [ "$DIM" == "2" ]; then
	case "$DIST" in
	"g")	RANDSOURCE="randn"
		PLAMBDA_SECOND="" ;;
	"u")	RANDSOURCE="randn"
		PLAMBDA_SECOND="2 sqrt / erf 3 sqrt *" ;;
	"c")	RANDSOURCE="randc"
		PLAMBDA_SECOND="0.01 *" ;;
	"l")	RANDSOURCE="randn randn randn randn 4 njoin"
		PLAMBDA_SECOND="x[0] x[1] * x[2] x[3] * - x[4] x[5] * x[6] x[7] * - join 2 sqrt /" ;;
	*)	echo -e "ERROR($0): unrecognized dist '$DIST'" >&2 false ;;
	esac
	PLAMBDA_FIRST="$RANDSOURCE $RANDSOURCE join pi sqrt * 2 * $GRAINSIZE *"
	PLAMBDA_SECOND="$PLAMBDA_SECOND $STRENGTH *"
else
	echo -e "ERROR($0): dimension \"$DIM\" not supported" >&2
	false
fi


echo -e "PLAMBDA_FIRST=\"$PLAMBDA_FIRST\"" >&2
echo -e "PLAMBDA_SECOND=\"$PLAMBDA_SECOND\"" >&2

plambda zero:${WIDTH}x${HEIGHT} "$PLAMBDA_FIRST" | blur $GRAINTYPE $GRAINSIZE | plambda - "$PLAMBDA_SECOND"
