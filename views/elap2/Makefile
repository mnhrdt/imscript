
C99 ?= clang

WFLAGS ?= -pedantic -Wall -Wextra -Wshadow

CFLAGS = $(WFLAGS) -g
CFLAGS = $(WFLAGS) -O3 -DNDEBUG

all: elap elap_order elap_init

elap: elap.c iio.o
	$(C99) $(CFLAGS) $< -o $@ iio.o -ljpeg -ltiff -lpng

elap_order: elap_order.c iio.o distance.o
	$(C99) $(CFLAGS) $< -o $@ iio.o distance.o -ljpeg -ltiff -lpng

elap_init: elap_init.c iio.o
	$(C99) $(CFLAGS) $< -o $@ iio.o -ljpeg -ltiff -lpng

iio.o: iio.c
	$(C99) $(CFLAGS) -c $< -o $@

distance.o: distance.c
	$(C99) $(CFLAGS) -DOMIT_DISTANCE_MAIN -c $< -o $@

test: elap elap_order elap_init
	./elap 0.49 300 lena.png mask.png result.png
	./elap_order 0.49 5 lena.png mask.png result_lex.png lex
	./elap_order 0.49 5 lena.png mask.png result_rand.png rand
	./elap_order 0.49 5 lena.png mask.png result_peel.png peel
	./elap_order 0.49 5 lena.png mask.png result_center.png center
	./elap_init 0.49 15 lena.png mask.png result_i0.png 0
	./elap_init 0.49 15 lena.png mask.png result_i127.png 127
	./elap_init 0.49 15 lena.png mask.png result_i255.png 255
	./elap_init 0.49 15 lena.png mask.png result_irand.png nan
	./elap_init 0.49 15 lena.png mask.png result_il.png lowres.png

clean:
	rm -rf iio.o elap elap_order elap_init result*.png distance.o
