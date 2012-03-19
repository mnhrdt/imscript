# compiler specific part (may be removed with minor damage)
#
ENABLE_GSL = yes
WFLAGS = -pedantic -Wall -Wextra -Wshadow -Wno-unused -Wno-array-bounds

CFLAGS = $(WFLAGS) -O3 -DNDEBUG
CFLAGS = $(WFLAGS)
CFLAGS = -g

#end of compiler specific part



SRCDIR = src
BINDIR = bin

SRCIIO = fftshift sterint plambda viewflow imprintf ntiply backflow unalpha imdim downsa flowarrows flowdiv fnorm imgstats qauto qeasy lrcat lk hs rgbcube iminfo setdim synflow vecstack ofc component faxpb faxpby iion flowgrad frakes_monaco_smith fillcorners colorflow lic deframe crosses
SRCFFT = gblur fft dct
ifeq ($(ENABLE_GSL), yes)
	SRCGSL = paraflow minimize
endif

IIOFLAGS = -ljpeg -ltiff -lpng
FFTFLAGS = -lfftw3f
GSLFLAGS = -lgsl -lgslcblas

# compiler detection hacks
# (because some compilers do not use the standard by default)
# TODO: move this stuff to a separate "hacks" file
ifeq ($(CC), cc)
	CC += -std=c99
endif
ifeq ($(CC), gcc)
	CC += -std=c99
endif
ifeq ($(CC), icc)
	CC += -std=c99
endif


# OS detection hacks
# TODO: move this stuff to a separate "portability" file
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CFLAGSIIO = $(CFLAGS) -D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=700
endif
ifeq ($(UNAME), Darwin)
	MPE := $(shell if test -d /opt/macports ; then echo yes ; fi)
	ifeq ($(MPE), yes)
		export C_INCLUDE_PATH := /opt/macports/include:$(C_INCLUDE_PATH)
		export LIBRARY_PATH := /opt/macports/lib:$(LIBRARY_PATH)
	endif
	MPE := $(shell if test -d /usr/X11 ; then echo yes ; fi)
	ifeq ($(MPE), yes)
		export C_INCLUDE_PATH := /usr/X11/include:$(C_INCLUDE_PATH)
		export LIBRARY_PATH := /usr/X11/lib:$(LIBRARY_PATH)
	endif
endif



SRC = $(SRCIIO) $(SRCFFT) $(SRCGSL)
PROGRAMS = $(addprefix $(BINDIR)/,$(SRC) flow_ms rgfield rgfields rgfieldst)


.PHONY: default
default: $(PROGRAMS)


$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS)

$(addprefix $(BINDIR)/,$(SRCFFT)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(addprefix $(BINDIR)/,$(SRCGSL)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(GSLFLAGS)

$(SRCDIR)/iio.o : $(SRCDIR)/iio.c $(SRCDIR)/iio.h
	$(CC) $(CFLAGSIIO) $(OFLAGS) -c $< -o $@

$(SRCDIR)/hs.o: $(SRCDIR)/hs.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_MAIN -c $< -o $@

$(SRCDIR)/gblur.o: $(SRCDIR)/gblur.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_GBLUR_MAIN -c $< -o $@

$(BINDIR)/flow_ms: $(addprefix $(SRCDIR)/,flow_ms.c gblur.o hs.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) -DUSE_MAINPHS $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(BINDIR)/rgfield: $(addprefix $(SRCDIR)/,rgfield.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(BINDIR)/rgfields: $(addprefix $(SRCDIR)/,rgfields.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(BINDIR)/rgfieldst: $(addprefix $(SRCDIR)/,rgfieldst.c gblur.o iio.o)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)


.PHONY: clean
clean:
	@rm -f $(PROGRAMS) $(SRCDIR)/*.o

.PHONY: zipdate
zipdate: clean
	(cd ..;tar --exclude-vcs -zchf imscript_`date +%Y_%m_%d_%H_%M`.tar.gz imscript)

.PHONY: zip
zip: clean
	(cd ..;tar --exclude-vcs -zchf imscript.tar.gz imscript)
