#!/bin/sh

PROGRAMS="blur fft iion imprintf plambda qeasy ransac siftu srmatch synflow pview vecov ifft lowe_join"

RM="rm -f"

$RM bin/iio.o
for i in $PROGRAMS; do
	$RM bin/$i
done
