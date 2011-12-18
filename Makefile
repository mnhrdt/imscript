CFLAGS += -std=c99

SRCDIR = src
BINDIR = bin

SRCIIO = plambda viewflow imprintf ntiply backflow unalpha imdim downsa flowarrows flowdiv fnorm imgstats qauto qeasy lrcat lk hs rgbcube iminfo
SRCFFT = gblur
SRCGSL = paraflow

SRC = $(SRCIIO) $(SRCFFT) $(SRCGSL)
PROGRAMS = $(addprefix $(BINDIR)/,$(SRC))
IIOFLAGS = $(SRCDIR)/iio.o -ljpeg -ltiff -lpng
FFTFLAGS = -lfftw3f
GSLFLAGS = -lgsl -lgslcblas


# OS detection hacks
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CFLAGS += -D_POSIX_C_SOURCE=200809L
endif
ifeq ($(UNAME), Darwin)
	MPE := $(shell if test -d /opt/macports ; then echo macport ; fi)
	ifeq ($(MPE), macport)
		export C_INCLUDE_PATH = /opt/macports/include
		export LIBRARY_PATH = /opt/macports/lib
	endif
endif

default: $(PROGRAMS)


$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@ $(IIOFLAGS)

$(addprefix $(BINDIR)/,$(SRCFFT)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@ $(IIOFLAGS) $(FFTFLAGS)

$(addprefix $(BINDIR)/,$(SRCGSL)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@ $(IIOFLAGS) $(GSLFLAGS)

$(SRCDIR)/iio.o : $(SRCDIR)/iio.c $(SRCDIR)/iio.h
	$(CC) $(CFLAGS) $(OFLAGS) -c $< -o $@

$(SRCDIR)/hs.o: $(SRCDIR)/hs.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_MAIN -c $< -o $@

$(SRCDIR)/gblur.o: $(SRCDIR)/gblur.c
	$(CC) $(CFLAGS) $(OFLAGS) -DOMIT_GBLUR_MAIN -c $< -o $@

$(BINDIR)/flow_ms: $(SRCDIR)/flow_ms.c $(BINDIR)/gblur.o $(BINDIR)/hs.o
	$(CC) $(CFLAGS) $(OFLAGS) -DUSE_MAINPHS $^ -o $@ $(IIOFLAGS) $(FFTFLAGS)


clean:
	rm -f $(PROGRAMS) $(BINDIR)/*.o

