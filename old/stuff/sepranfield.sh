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

case "$DIST" in
"g")	RANDSOURCE="randn"
	PLAMBDA_SECOND="" ;;
"u")	RANDSOURCE="randn"
	PLAMBDA_SECOND="2 sqrt / erf 3 sqrt *" ;;
"c")	RANDSOURCE="randc"
	PLAMBDA_SECOND="0.01 *" ;;
"l")	RANDSOURCE="randn randn randn randn 4 njoin"
	PLAMBDA_SECOND="4 nsplit * >1 * <1 - 2 sqrt /" ;;
*)	echo -e "ERROR($0): unrecognized dist '$DIST'" >&2 false ;;
esac

RANDSOURCES=`yes $RANDSOURCE|sed ${DIM}q|tr '\n' ' '`
PLAMBDA_FIRST="$RANDSOURCES $DIM njoin pi sqrt * 2 * $GRAINSIZE *"
PLAMBDA_SECOND="$PLAMBDA_SECOND $STRENGTH *"


#echo -e "PLAMBDA_FIRST=\"$PLAMBDA_FIRST\"" >&2
#echo -e "PLAMBDA_SECOND=\"$PLAMBDA_SECOND\"" >&2

plambda zero:${WIDTH}x${HEIGHT} "$PLAMBDA_FIRST" | blur $GRAINTYPE $GRAINSIZE | plambda - "$PLAMBDA_SECOND"
