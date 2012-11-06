#!/bin/bash

# usage:
# graydsift.sh frame.pgm  > frame.sift.txt

SIFT_EXECUTABLE_LOWE=sift_executable_lowe
SIFT_EXECUTABLE_HESS=sift_executable_hess
SIFT_EXECUTABLE_CAO=lowe_sift_yu_v1
SIFT_EXECUTABLE_ZEROFROG=sift_executable_zerofrog
SIFT_EXECUTABLE_VEDALDI=sift_executable_vedaldi
SIFT_EXECUTABLE_SURF=sift_executable_surf
SIFT_EXECUTABLE_RIDA=sift_executable_rida

set -e

IFILE=$2
IMETH=$1

usgexit() {
	echo -e "usage:\n\t `basename $0` [lowe,cao,zero,surf] frame.pgm > frame.sift" >&2
	exit 1
}

testex() {
	if which $1 > /dev/null ; then
		echo > /dev/null
	else
		echo "ERROR: executable file $1 not available" >&2
		exit 1
	fi
}

sift_from_cao() {
	testex $SIFT_EXECUTABLE_CAO
	testex flprintasc
	TPD=`mktemp -d /tmp/dsiftie.XXXXXX` || exit 1
	$SIFT_EXECUTABLE_CAO $1 $TPD/sift
	flprintasc $TPD/sift
	rm -rf $TPD
}

sift_from_lowe() {
	testex $SIFT_EXECUTABLE_LOWE
	testex lowe_join
	testex plambda
	testex siftu
	qeasy 0 255 $1 | $SIFT_EXECUTABLE_LOWE | lowe_join | siftu aff 0 1 0 1 0 0 | sort -n -k 3
}

sift_from_zerofrog() {
	testex $SIFT_EXECUTABLE_ZEROFROG
	testex lowe_join
	testex plambda
	testex siftu
	qeasy 0 255 $1 | $SIFT_EXECUTABLE_ZEROFROG | lowe_join | siftu aff 0 1 0 1 0 0 | sort -n -k 3
}

sift_from_vedaldi() {
	testex $SIFT_EXECUTABLE_VEDALDI
	TPD=`mktemp -d /tmp/dsiftie.XXXXXX` || exit 1
	$SIFT_EXECUTABLE_VEDALDI -o $TPD/s $1
	cat $TPD/s
	rm -rf $TPD
}

sift_from_hess() {
	testex $SIFT_EXECUTABLE_HESS
	TPD=`mktemp -d /tmp/dsiftie.XXXXXX` || exit 1
	$SIFT_EXECUTABLE_HESS -x -o $TPD/s $1
	cat $TPD/s
	rm -rf $TPD
}

sift_from_rida() {
	testex $SIFT_EXECUTABLE_RIDA
	TPD=`mktemp -d /tmp/dsiftie.XXXXXX` || exit 1
	mkdir -p $TPD/images/tmp
	cp $1 $TPD
	cd $TPD
	$SIFT_EXECUTABLE_RIDA -s --post --ck 100 --si `basename $1` --so x.ssift 1>&2
	cd - >/dev/null
	tail -n +2 $TPD/x.ssift
	rm -rf $TPD
}


sift_from_surf() {
	testex $SIFT_EXECUTABLE_SURF
	testex qnm
	TPD=`mktemp -d /tmp/dsiftie.XXXXXX` || exit 1
	qnm pgmbin < $2 > $TPD/g
	case "$1" in
	"normal")   SURFO="" ;;
	"128")      SURFO="-e" ;;
	"upright")  SURFO="-u" ;;
	"u128")     SURFO="-u -e" ;;
	"double")   SURFO="-d" ;;
	"d128")     SURFO="-d -e" ;;
	"dupright") SURFO="-d -u" ;;
	"du128")    SURFO="-d -u -e" ;;
	*) usgexit;
	esac
	$SIFT_EXECUTABLE_SURF -i $TPD/g -o $TPD/s $2 -q $SURFO
	tail -n +3 $TPD/s
	rm -rf $TPD
}


# check input
if [ $# != "2" ]; then
	usgexit
fi

if [ "`file -Lb $IFILE|cut -c1-10`" != "Netpbm PGM" ]; then
	echo "input file \"$IFILE\" is not in PGM format" >&2
	exit 1
fi


case "$IMETH" in
"lowe")
	PRUN_SIFT=sift_from_lowe
	;;
"hess")
	PRUN_SIFT=sift_from_hess
	;;
"rida")
	PRUN_SIFT=sift_from_rida
	;;
"cao")
	PRUN_SIFT=sift_from_cao
	;;
"zero")
	PRUN_SIFT=sift_from_zerofrog
	;;
"vedaldi")
	PRUN_SIFT=sift_from_vedaldi
	;;
"surf")
	PRUN_SIFT="sift_from_surf normal"
	;;
"surf128")
	PRUN_SIFT="sift_from_surf 128"
	;;
"surfu")
	PRUN_SIFT="sift_from_surf upright"
	;;
"surfu128")
	PRUN_SIFT="sift_from_surf u128"
	;;
"surfd")
	PRUN_SIFT="sift_from_surf double"
	;;
"surfd128")
	PRUN_SIFT="sift_from_surf d128"
	;;
"surfdu")
	PRUN_SIFT="sift_from_surf dupright"
	;;
"surfdu128")
	PRUN_SIFT="sift_from_surf du128"
	;;
*)
	usgexit;
esac

$PRUN_SIFT $IFILE
