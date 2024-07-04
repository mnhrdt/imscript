CFLAGS ?= -O3
LDLIBS += -lm -lfftw3f

OBJ = src/iio.o src/fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa autotrim iion mediator     \
      redim colormatch eucdist nonmaxsup gntiply idump warp heatd imhalve    \
      ppsmooth mdither mdither2 rpctk getbands pixdump bandslice points      \
      columnize autocorr
      #geomedian carve

BIN := $(addprefix bin/,$(BIN))

default: $(BIN) bin/cpu_term bin/rpcflip_term

bin/%  : src/%.o $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)


# CONFIGURABLE DEPENDENCIES
# XXX: comment (or not) the following lines to disable (enable) image formats
# Note: just comment the lines, do NOT change the ones to zeros!

ENABLE_PNG  = 1
ENABLE_TIFF = 1
ENABLE_JPEG = 1
ENABLE_WEBP = 1
ENABLE_HEIF = 1
ENABLE_PGSL = 1

# CAVEAT: if you want to use HDF5, make sure that no "mpich" packages
# are installed on your computer.  If they are, all programs that link
# to libdf5 become really slow due to dynamic linking initialization.
#
ENABLE_HDF5 = 1





#
# END OF CORE PART, THE REST OF THIS FILE IS NOT ESSENTIAL
#


# unit test (end-to-end)
ifeq ($(MAKECMDGOALS),test)
LDLIBS := $(filter-out -lfftw3f,$(LDLIBS))
endif
test: bin/plambda bin/imprintf
	echo $(MAKECMDGOALS)
	bin/plambda zero:512x512 "randg randg randl randg randg 5 njoin" \
	| bin/plambda "split rot del hypot" | bin/imprintf "%s%e%r%i%a\n" \
	| grep -q 3775241.440161.730120.0037888911.8794


# build tutorial
tutorial: default
	cd doc/tutorial/f && env PATH=../../../bin:$(PATH) $(SHELL) build.sh

# build manpages
manpages: default
	cd doc && $(SHELL) rebuild_manpages.sh

# deploy stagit
stagit:       ; cd doc/stagit && sh build.sh
stagit-clean: ; cd doc/stagit && sh clean.sh



# exotic targets, not built by default
# FTR: interactive tools, require X11 or freeglut
# MSC: "misc" tools, or those requiring GSL

#OBJ_FTR = src/iio.o src/ftr/ftr.o src/ftr/egm96.o src/ftr/fancy_image.o src/ftr/ccpu.o
OBJ_FTR = src/iio.o src/ftr/ftr.o src/ftr/fancy_image.o src/ftr/ccpu.o

BIN_FTR = $(shell test -f src/ftr/TARGETS && cat src/ftr/TARGETS)
BIN_MSC = $(shell test -f src/misc/TARGETS && cat src/misc/TARGETS)

BIN_FTR := $(addprefix bin/,$(BIN_FTR))
BIN_MSC := $(addprefix bin/,$(BIN_MSC))

LDLIBS_FTR = $(LDLIBS)
LDLIBS_MSC = $(LDLIBS) -lgsl -lgslcblas

OBJ_ALL = $(OBJ) $(OBJ_FTR)
BIN_ALL = $(BIN) $(BIN_FTR) $(BIN_MSC)

full  : default ftr misc
ftr   : $(BIN_FTR)
misc  : $(BIN_MSC)


bin/% : src/ftr/%.o $(OBJ_FTR)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS_FTR)

bin/% : src/misc/%.o $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS_MSC)


ifdef FTR_GLUT
LDLIBS_FTR += -lGL -lglut
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBTIFF
src/ftr/ftr.o : CFLAGS += -DFTR_BACKEND=\'f\'
else
LDLIBS_FTR += -lX11
endif



# this is the end of the "regular" targets.  What follows are fancy build
# targets, or bizarre graphics/terminal interfaces.

# static libs need all to be explicitly pulled (no recursion)
# (NOTE: assumes you want jpeg, tiff, png, webp support)
STALIBS = -lm -lfftw3f -ltiff -ljpeg -lpng -lwebp -ljbig -lz -llzma -ldeflate -lzstd -lX11 -lxcb -lm -lXau -lXdmcp $(STALIBSX)
STALIBSX = -lX11 -lxcb -lXau -lXdmcp
# uncomment the following two lines for static compilation
#LDFLAGS = -static
#LDLIBS = $(STALIBS)


# single, fat, busybox-like executable
BINOBJ = $(BIN:bin/%=src/%.o) $(BIN_FTR:bin/%=src/ftr/%.o) src/ftr/ftr.o
#L = -lfftw3f -lpng -ltiff -ljpeg -llzma -lz -lm -ljbig $(LDLIBS) -lm -lpthread -lgsl -lgslcblas -lzstd
L = -lfftw3f -lpng -ltiff -ljpeg -llzma -lz -ldeflate -lm -ljbig -lwebp -lpthread -lzstd -lX11 -lxcb -ldl -lgsl -lgslcblas -lhdf5_serial -ldl -lm -lXau -lpthread -lXdmcp -lzstd -lz -lsz -laec
bin/im.static : src/im.o $(BINOBJ) $(OBJ) src/misc/overflow.o
	$(CC) $(LDFLAGS) -static -Wl,--allow-multiple-definition -o $@ $^ $L

BINOBJ2 = $(BIN:bin/%=src/%.o) $(BIN_FTR:bin/%=src/ftr/%.o) src/ftr/ftr.o
L2 = $(LDLIBS) $(LDLIBS_FTR) -lgsl
bin/im : src/im.o $(BINOBJ2) $(OBJ) src/misc/overflow.o
	$(CC) $(LDFLAGS) -Wl,--allow-multiple-definition -o $@ $^ $(L2)


# some ftr executables, but compiled for the terminal backend
OBJ_FTR_TERM = src/ftr/ftr_term.o $(filter-out src/ftr/ftr.o,$(OBJ_FTR))
bin/%_term : src/ftr/%.o $(OBJ_FTR_TERM)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

OBJ_FTR_SIX = src/ftr/ftr_sixel.o $(filter-out src/ftr/ftr.o,$(OBJ_FTR))
bin/%_six : src/ftr/%.o $(OBJ_FTR_SIX)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)







# bureaucracy
clean: ; @$(RM) $(BIN_ALL) bin/im src/*.o src/ftr/*.o src/misc/*.o
.PHONY: default full ftr misc clean tutorial manpages
.PRECIOUS: %.o


# hack (compatibility flags for old compilers)
#
# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
#ifeq (,$(shell $(CC) $(CFLAGS) -dM -E -< /dev/null | grep __STDC_VERSION_))
#CFLAGS := $(CFLAGS) -std=gnu99
#endif



#
# setup variables for configurable dependencies (nothing interesting here)
#

ifdef ENABLE_PNG
LDLIBS += -lpng
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBPNG
endif

ifdef ENABLE_JPEG
LDLIBS += -ljpeg
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBJPEG
endif

ifdef ENABLE_WEBP
LDLIBS += -lwebp
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBWEBP
endif

ifdef ENABLE_HEIF
LDLIBS += -lheif
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBHEIF
endif

ifdef ENABLE_TIFF
LDLIBS += -ltiff
src/iio.o: CPPFLAGS += -DI_CAN_HAS_LIBTIFF
else
src/fancy_image.o: CPPFLAGS += -DFANCY_IMAGE_DISABLE_TIFF
# note that disabling tiff kills the whole fancy_image stuff
endif

ifdef ENABLE_HDF5
LDLIBS += $(shell pkg-config hdf5 --libs --silence-errors || echo -lhdf5)
src/iio.o: CPPFLAGS+= -DI_CAN_HAS_LIBHDF5 `pkg-config hdf5 --cflags 2>/dev/null`
# yes, the hdf5 compile-time configuration is a bit fucked up, but this goes
# well with the rest of that library.
endif

ifdef ENABLE_PGSL
bin/plambda: LDLIBS += -lgsl
src/plambda.o: CPPFLAGS += -DPLAMBDA_WITH_GSL
bin/mathieu: LDLIBS += -lgsl
endif




# create and include non-standard file dependences
# (this is needed because some .c files include other .c files directly)
# it's not a big deal, you can comment this lines if they fail
D = src src/ftr src/misc
src/dep.mk : ; for i in $D;do $(CC) -MM $$i/*.c|sed "\:^[^ ]:s:^:$$i/:g";done>$@
-include src/dep.mk
