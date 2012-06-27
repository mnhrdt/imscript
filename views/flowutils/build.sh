CC="cc -std=c99"
LFLAGS="-lpng -ltiff"

$CC -D_XOPEN_SOURCE=700 -c iio.c
$CC backflow.c iio.o -o backflow $LFLAGS
$CC flowarrows.c iio.o -o flowarrows $LFLAGS
$CC flowinv.c iio.o -o flowinv $LFLAGS
$CC viewflow.c iio.o -o viewflow $LFLAGS
$CC qeasy.c iio.o -o qeasy $LFLAGS
