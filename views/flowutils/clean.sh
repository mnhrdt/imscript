PROGRAMS="backflow flowarrows flowinv viewflow qeasy plambda blur"

rm -f bin/iio.o
for i in $PROGRAMS; do
	rm -f bin/$i
done
rm -f test/{hom_arrows.png,hom_colors.png,hom_lena.png,ihom.flo,ihom_lena.png,ihom_arrows.png}
