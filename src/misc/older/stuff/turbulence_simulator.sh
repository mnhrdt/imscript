#!/bin/bash

set -e

export PATH=$PATH:/bin:/usr/bin

if [ "$#" != "4" ]; then
	echo -e "usage:\n\t$0 "\
	"in.png nframes grain strength"
	#1      2       3     4
	false
fi


INFILE=$1
NFRAMES=$2
GRAINSIZE=$3
STRENGTH=$4
WIDTH=`imprintf %w $INFILE`
HEIGHT=`imprintf %h $INFILE`

for i in `seq 1 $NFRAMES`; do
	sepranfield.sh $WIDTH $HEIGHT 2 g $GRAINSIZE g $STRENGTH $[42+i] | backflow - $INFILE out_$i.png &
done
wait
