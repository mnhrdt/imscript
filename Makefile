#CFLAGS ?= -march=native -O3 -DNDEBUG -Wall -Wno-unused
#CFLAGS ?= -march=native -O3 -Wall -Wextra -Wno-unused  -fsanitize=addre
#CFLAGS ?= -g -Wall -Wextra -Wno-unused #-fsanitize=address
CFLAGS ?= -march=native -O3 -DNDEBUG -Wall -Wno-unused
LDLIBS += -ljpeg -ltiff -lpng -lz -lfftw3f -lm #-lgdal

OBJ = src/iio.o src/fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa autotrim iion mediator     \
      redim colormatch eucdist nonmaxsup gntiply idump warp heatd imhalve    \
      ppsmooth mdither mdither2 rpctk

BIN := $(addprefix bin/,$(BIN))

default: $(BIN) bin/cpu_term bin/rpcflip_term

bin/%  : src/%.o $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)


# END OF CORE PART, THE REST OF THIS FILE IS NOT ESSENTIAL

# unit test (end-to-end)
ifeq ($(MAKECMDGOALS),test)
LDLIBS := $(filter-out -lfftw3f,$(LDLIBS))
endif
test: bin/plambda bin/imprintf
	echo $(MAKECMDGOALS)
	bin/plambda zero:512x512 "randg randg randl randg randg 5 njoin" \
	| bin/plambda - "split rot del hypot" | bin/imprintf "%s%e%r%i%a\n" \
	| grep -q 3775241.440161.730120.0037888911.8794


# hack (compatibility flags for old compilers)
#
# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
ifeq (,$(shell $(CC) $(CFLAGS) -dM -E -< /dev/null | grep __STDC_VERSION_))
CFLAGS := $(CFLAGS) -std=gnu99
endif



# exotic targets, not built by default
# FTR: interactive tools, require X11 or freeglut
# MSC: "misc" tools, or those requiring GSL

OBJ_FTR = src/iio.o src/ftr/ftr.o src/ftr/egm96.o src/ftr/fancy_image.o

BIN_FTR = $(shell cat src/ftr/TARGETS)
BIN_MSC = $(shell cat src/misc/TARGETS)

BIN_FTR := $(addprefix bin/,$(BIN_FTR))
BIN_MSC := $(addprefix bin/,$(BIN_MSC))

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

# single, fat, busybox-like executable
BINOBJ = $(BIN:bin/%=src/%.o) #$(BIN_FTR:bin/%=src/ftr/%.o)
L = -lfftw3f -lpng -ltiff -ljpeg -llzma -lz -lm -ljbig
bin/im : src/im.o $(BINOBJ) $(OBJ) src/misc/overflow.o
	$(CC) $(LDFLAGS) -static -Wl,--allow-multiple-definition -o $@ $^ $L

# some ftr executable, but compiled for the terminal backend
OBJ_FTR_TERM = src/ftr/ftr_term.o $(filter-out src/ftr/ftr.o,$(OBJ_FTR))
bin/%_term : src/ftr/%.o $(OBJ_FTR_TERM)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# bureaucracy
clean: ; @$(RM) $(BIN_ALL) bin/im src/*.o src/ftr/*.o src/misc/*.o
.PHONY: default full ftr misc clean
.PRECIOUS: %.o

# non-standard file dependences
# (this is needed because some .c files include other .c files directly)
DIRS = src src/ftr src/misc
.deps.mk: ; for i in $(DIRS);do cc -MM $$i/*.c|sed "\:^[^ ]:s:^:$$i/:g";done>$@
-include .deps.mk
