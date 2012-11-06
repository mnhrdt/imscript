#!/bin/sh


CC="c99"

#other possibilities for the compiler:
# CC="clang"
# CC="sunc99"
# CC="icc -std=c99"
# CC="gcc -std=c99"
# CC="gcc -std=c99 -pedantic -Wall -Wextra -Werror -Wno-unused"


PROGRAMS="blur fft iion imprintf plambda qeasy ransac siftu srmatch synflow pview vecov lowe_join"


LFLAGS="-lpng -ltiff -lfftw3f"


$CC -D_XOPEN_SOURCE=700 -c src/iio.c -o bin/iio.o
for i in $PROGRAMS; do
	$CC src/$i.c bin/iio.o -o bin/$i $LFLAGS
done
ln -sf fft bin/ifft
