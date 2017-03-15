# compiler specific part (may be removed with minor damage)
#
ENABLE_GSL = yes
WFLAGS=
WFLAGS = -pedantic -Wall -Wextra -Wshadow -Wstrict-prototypes -Wno-unused -Wno-parentheses -Wno-sign-compare -Werror -Wno-error=format -Wno-error=overflow
WFLAGS = -pedantic -Wall -Wextra -Wshadow -Wstrict-prototypes -Wno-unused -Wno-sign-compare
WFLAGS=

OFLAGS = $(WFLAGS)
OFLAGS = -march=native -O3

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) $(CFLAGS) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif



SRCIIO = fftshift plambda sterint viewflow imprintf ntiply backflow unalpha imdim downsa flowarrows flowdiv fnorm imgstats qauto qeasy lrcat tbcat lk hs rgbcube iminfo setdim synflow vecstack ofc component faxpb faxpby iion flowgrad fillcorners colorflow lic deframe crosses crop angleplot closeup hrezoom upsa veco vecov vecov_lm flowinv ghisto shuntingyard rpc overpoints periodize rpcflow ransac genk cgi zeropad siftu pview homfilt rpchfilt rpc_errfilt uncrop maptp rpcparcheck rpc_errsingle rpc_errpair cline rpc_angpair rpc_curvpair chisto fftper srmatch croparound zoombil flowh harris lgblur rpc_eval cutrecombine cdeint imgerr amle elap imspread replicate disp_to_corr printmask colormesh colormeshh ijmesh elap3 sphereheights fmsrA elap_rec amle_rec palette radphar aronsson43 scheme_plap flowjac elap_recsep aronsson11 fill_bill fill_rect rpc_epicyl starfield morsi gharrows ipol_watermark fontu fontu2 cglap flownop pairsinp pairhom poisson_rec cgpois cgpois_rec isoricci lapbediag lapcolo simplest_inpainting lapbediag_sep cldmask plyflatten metatiler hview dither ditheru histeq8 thinpa_recsep really_simplest_inpainting bmms perms censust satproj mnehs mnehs_ms aff3d amle_recsep elevate_matches elevate_matcheshh pmba pmba2 frustumize ghough ghough2 tdip ihough2 alphadots posmax ppsboundary homdots awgn strt remove_small_cc distance tregistration tcregistration registration tvint shadowcast raddots simpois nnint bdint contihist graysing unshadow sfblur overflow ppsmooth inppairs mima utm autholes plyroads startracker homwarp tiffu fancy_evals rpc_warpab rpc_warpabt rpc_mnehs rpc_pm rpc_pmn bandhisto polygonify dveco dvecov_lm drawroads bayerparts bayerdots vecoh
SRCFFT = gblur fft dct blur lgblur2 lure lgblur3 testgblur lures dht
ifeq ($(ENABLE_GSL), yes)
	SRCGSL = paraflow minimize
endif

# swap next two lines for libraw support
IIOFLAGS = -ljpeg -ltiff -lpng -lm -lraw -lHalf -lIex -lIlmImf -fopenmp
IIOFLAGS = -ljpeg -ltiff -lpng -lm
FFTFLAGS = -lfftw3f
GSLFLAGS = -lgsl -lgslcblas

# swap next two lines for libraw support
IIO=../iio/contrib/cmake/libIIOLIB.a
IIO=src/iio.o


SRC = $(SRCIIO) $(SRCFFT) $(SRCGSL)
PROGRAMS = $(addprefix bin/,$(SRC) flow_ms rgfield rgfields rgfieldst elap2 flambda fancy_zoomout fancy_downsa fancy_crop)


.PHONY: default
default: $(PROGRAMS)

$(addprefix bin/,depend) : src/

$(addprefix bin/,$(SRCIIO)) : bin/% : src/%.c $(IIO)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

$(addprefix bin/,$(SRCFFT)) : bin/% : src/%.c $(IIO)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(addprefix bin/,$(SRCGSL)) : bin/% : src/%.c $(IIO)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(GSLFLAGS)

src/iio.o : src/iio.c src/iio.h
	$(CC) $(CFLAGS) $(OFLAGS) -c $< -o $@

src/hs.o: src/hs.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_MAIN -c $< -o $@

src/lk.o: src/lk.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_MAIN -c $< -o $@

src/gblur.o: src/gblur.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_GBLUR_MAIN -c $< -o $@

src/flow_ms.o: src/flow_ms.c
	$(CC) $(CFLAGS) $(OFLAGS) -c $< -o $@

src/flowarrows.o: src/flowarrows.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_MAIN -c $< -o $@

#ifndef OMIT_BLUR_MAIN
#define MAIN_BLUR
#endif
src/distance.o: src/distance.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_DISTANCE_MAIN -c $< -o $@

bin/flow_ms: $(addprefix src/,flow_ms.c gblur.o hs.o lk.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) -DUSE_MAINAPI $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

bin/rgfield: $(addprefix src/,rgfield.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

bin/rgfields: $(addprefix src/,rgfields.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

bin/rgfieldst: $(addprefix src/,rgfieldst.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

bin/elap2: $(addprefix src/,elap2.c distance.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

bin/flambda: $(addprefix src/,flambda.c fancy_image.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

bin/fancy_zoomout: $(addprefix src/,fancy_zoomout.c fancy_image.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

bin/fancy_downsa: $(addprefix src/,fancy_downsa.c fancy_image.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

bin/fancy_crop: $(addprefix src/,fancy_crop.c fancy_image.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

.PHONY: clean
clean:
	@rm -f $(PROGRAMS) src/*.o

.PHONY: zipdate
zipdate: clean
	(cd ..;tar --exclude-vcs -zchf imscript_`date +%Y_%m_%d_%H_%M`.tar.gz imscript)

.PHONY: zip
zip: clean
	(cd ..;tar --exclude-vcs -zchf imscript.tar.gz imscript)

.PHONY: depend
depend:
	$(CC) -MM src/*.c 2>/dev/null|sed '/^[^ ]/s/^/src\//'|sed 's/\.o:/:/'>src/dependencies

-include src/dependencies


# TODO: implement the following behavior without the touch or other ugly hacks
#
#bin/plambda: src/plambda.c src/iio.o
#	$(CC) $(CFLAGS) -o $@ $^ $(IIOFLAGS)
#bin/plambda: src/plambda.o src/iio.o
#	$(CC) -o $@ $^ $(IIOFLAGS)
#
#src/plambda.o: src/random.c
#	@touch src/plambda.c




test: bin/plambda bin/imprintf
	bin/plambda zero:512x512 "rand rand rand rand rand 5 njoin" | bin/plambda - "split del 0 >" | bin/imprintf "%s\n" | grep -q 131148
