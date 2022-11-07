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
#ENABLE_WEBP = 1
#ENABLE_HEIF = 1
#ENABLE_PGSL = 1

# NOTE **WARNING**
#
# You can uncomment the line below to enable support for HDF5 files.
# However, notice that libHDF5 is an extremely user-hostile library and
# linking to it has major drawbacks in exchange of being able to open the
# stupid hdf5 files.  The real nightmare starts after version 1.10.3 of that
# library, which is based on the "fabric" framework.  Thus, linking against
# hdf5 turns each program into an "hdf5 application" that will perform
# several "checks" each time the program is run.  These checks are triggered
# by the mere dynamic loading of the library and cannot be disabled (even if
# the program does not use the library at all!).  On my laptop, this takes
# TWO HUNDRED FUCKING MILLISECONDS for each run of each program.  This is
# easily way longer than most actual image-processing computations, slowing
# all imscript-based pipelines to a crawl.  It's not that these checks do
# anything useful, either : if you strace the runtime linking of hdf5, you'll
# see that most of this time is spent in 961 calls to clock_nanosleep(2).
#
# My disappointment is immeasurable, and my day is ruined.
#
# On future versions of imscript, HDF5 support will be enabled by means of a
# middleman library that dynamically links libhdf5 upon demand.  Thus, only
# the first reading of an actual hdf5 image will be slowed-down.  But this
# is a lot of work and it will have to wait.
#
#ENABLE_HDF5 = 1





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



# exotic targets, not built by default
# FTR: interactive tools, require X11 or freeglut
# MSC: "misc" tools, or those requiring GSL

OBJ_FTR = src/iio.o src/ftr/ftr.o src/ftr/egm96.o src/ftr/fancy_image.o

BIN_FTR = $(shell cat src/ftr/TARGETS)
BIN_MSC = $(shell cat src/misc/TARGETS)

BIN_FTR := $(addprefix bin/,$(BIN_FTR))
BIN_MSC := $(addprefix bin/,$(BIN_MSC))

LDLIBS_FTR = $(LDLIBS) -lGL -lglut
LDLIBS_FTR = $(LDLIBS) -lX11
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


# this is the end of the "regular" targets.  What follows are fancy build
# targets, or bizarre graphics/terminal interfaces.


# single, fat, busybox-like executable
BINOBJ = $(BIN:bin/%=src/%.o) $(BIN_FTR:bin/%=src/ftr/%.o) src/ftr/ftr.o
#L = -lfftw3f -lpng -ltiff -ljpeg -llzma -lz -lm -ljbig $(LDLIBS) -lm -lpthread -lgsl -lgslcblas -lzstd
L = -lfftw3f -lpng -ltiff -ljpeg -llzma -lz -lm -ljbig -lwebp -lpthread -lzstd -lX11 -lxcb -ldl -lgsl -lgslcblas -lhdf5_serial -ldl -lm -lXau -lpthread -lXdmcp -lzstd -lz -lsz -laec
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
endif




# create and include non-standard file dependences
# (this is needed because some .c files include other .c files directly)
# it's not a big deal, you can comment this lines if they fail
DIRS = src src/ftr src/misc
.deps.mk:;for i in $(DIRS);do $(CC) -MM $$i/*.c|sed "\:^[^ ]:s:^:$$i/:g";done>$@
-include .deps.mk
