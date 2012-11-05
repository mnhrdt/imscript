CC="c99"
LFLAGS="-lpng -ltiff -lfftw3f"

PROGRAMS="blur fft iion imprintf plambda qeasy ransac siftu srmatch synflow"

$CC -D_XOPEN_SOURCE=700 -c src/iio.c -o bin/iio.o
for i in $PROGRAMS; do
	$CC src/$i.c bin/iio.o -o bin/$i $LFLAGS
done
