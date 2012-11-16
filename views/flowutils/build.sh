CC="cc -std=c99"
LFLAGS="-lpng -ltiff -lfftw3f"

PROGRAMS="backflow flowarrows flowinv viewflow qeasy plambda blur"

$CC -D_XOPEN_SOURCE=700 -c src/iio.c -o bin/iio.o
for i in $PROGRAMS; do
	$CC src/$i.c bin/iio.o -o bin/$i $LFLAGS
done
