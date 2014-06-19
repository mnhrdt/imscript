#!/bin/bash

# size of each component
GEOMETRY="3712 3202"

# SIFT parameters: sift threshold and sift matching distance
SSEUIL=290
SDIST=20

# tiling parameters: tiles accross and down
TILESX=11
TILESY=7

# ransac parameters
RNTRIALS=1000
RDIST=0.9
RMININ=30

# vector field regularity parameter
DELTA=100




# abbreviations
HIDENTITY="1 0 0 0 1 0 0 0 1"


set -e

usgexit() {
	echo -e "usage:\n\t$0 negative.tiff" >&2
	false
}

if [ "$#" != "1" ]; then
	usgexit
fi


FILE=$1
BASENAME=`basename $FILE .tiff`
PREFIX=`echo $FILE|sed 's/.tiff$//'`
REFILE=$PREFIX.tiff

if [ "$REFILE" != "$FILE" ]; then
	usgexit
fi

echo "FILE=$FILE"
echo "OUTPUT_DIRECTORY=$PREFIX"

OD=$PREFIX

mkdir -p $OD
cd $OD

ln -sf ../$BASENAME.tiff X.tiff

test ! -f G2.tiff && \
	cut3 X.tiff G0.tiff G1.tiff G2.tiff

for i in 0 1 2; do
	test ! -f G$i.png && \
		quantize 0 65535 G$i.tiff G$i.png
done

ln -sf G1.png g1.png
for i in 0 2; do
	test ! -f g$i.png && \
		registration G1.png G$i.png g$i.png
done

for i in 0 1 2; do
	test ! -f g$i.sift && \
		unisift.sh zerofrog g$i.png >g$i.sift
done

for i in 0 2; do
	test ! -f pairs.1$i && \
		siftu pairg $SSEUIL g1.sift g$i.sift $SDIST $SDIST pairs.1$i
	test ! -f vpairs.1$i.png && \
		pview pairs $HIDENTITY $GEOMETRY < pairs.1$i > vpairs.1$i.png
done

test ! -f crop.0 && \
	tiles $GEOMETRY $TILESX $TILESY | while read i a b c d ; do
		echo $a $b $c $d > crop.$i
	done

NTILES=`ls crop.*|wc -l`

for d in 0 2; do
	for i in `seq 0 $[NTILES-1]`; do
		test ! -f pairs.1$d.$i && \
			paircrop `cat crop.$i` < pairs.1$d > pairs.1$d.$i
	done
done

for d in 0 2; do
	for i in `seq 0 $[NTILES-1]`; do
		test ! -f hom.1$d.$i && \
			ransac hom $RNTRIALS $RDIST $RMININ hom.1$d.$i < pairs.1$d.$i
	done
done

for d in 0 2; do
	test ! -f flow_1$d.d$DELTA.tiff && \
		for i in `seq 0 $[NTILES-1]`; do paste crop.$i hom.1$d.$i ; done |
			combine_homographies $GEOMETRY $DELTA flow_1$d.d$DELTA.tiff
done

for d in 0 2; do
	test ! -f rg$d.$DELTA.png && \
		backflow flow_1$d.d$DELTA.tiff g$d.png rg$d.d$DELTA.png
done

test ! -f resultat.d$DELTA.png && \
	join3 rg2.d$DELTA.png g1.png rg0.d$DELTA.png resultat.d$DELTA.png


echo "everything done PWD=`pwd`"

# vim:set ts=4 sw=4:
