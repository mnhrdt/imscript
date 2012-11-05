RM="rm -f"
PROGRAMS="blur fft iion imprintf plambda qeasy ransac siftu srmatch synflow"

$RM bin/iio.o
for i in $PROGRAMS; do
	$RM bin/$i
done
