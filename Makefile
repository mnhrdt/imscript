SRCDIR = src
BINDIR = bin

SRCIIO = plambda viewflow imprintf ntiply backflow unalpha imdim downsa flowarrows flowdiv fnorm imgstats qauto qeasy lrcat lk hs rgbcube iminfo setdim synflow vecstack ofc component faxpb faxpby iion flowgrad
SRCFFT = gblur fft dct
SRCGSL = paraflow minimize

IIOFLAGS = -ljpeg -ltiff -lpng
FFTFLAGS = -lfftw3f
GSLFLAGS = -lgsl -lgslcblas

# compiler detection hacks (some compilers do not use the standard by default)
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
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CFLAGSIIO = $(CFLAGS) -D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=700
endif
ifeq ($(UNAME), Darwin)
	MPE := $(shell if test -d /opt/macports ; then echo macport ; fi)
	ifeq ($(MPE), macport)
		export C_INCLUDE_PATH = /opt/macports/include
		export LIBRARY_PATH = /opt/macports/lib
	endif
endif



SRC = $(SRCIIO) $(SRCFFT) $(SRCGSL)
PROGRAMS = $(addprefix $(BINDIR)/,$(SRC) flow_ms)


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


.PHONY: clean
clean:
	@rm -f $(PROGRAMS) $(SRCDIR)/*.o

.PHONY: zipdate
zipdate: clean
	(cd ..;tar --exclude-vcs -zchf imscript_`date +%Y_%m_%d_%H_%M`.tar.gz imscript)

.PHONY: zip
zip: clean
	(cd ..;tar --exclude-vcs -zchf imscript.tar.gz imscript)
