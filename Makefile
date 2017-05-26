# user-editable configuration (with reasonable defaults)
CFLAGS ?= -march=native -Os

# required libraries
LDLIBS = -ljpeg -ltiff -lpng -lz -lfftw3f -lm

# variables
OBJ = src/iio.o src/fancy_image.o
LIB = src/libiio.a
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa iion                       \

BIN := $(addprefix bin/,$(BIN))

# default target
default: $(BIN)

# rules
bin/% : src/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(LIB) : $(LIB)($(OBJ))

# bureaucracy
clean:
	$(RM) $(OBJ) $(LIB) $(BIN)
.PHONY: default clean

# END OF CORE PART, THE REST OF THIS FILE IS NOT ESSENTIAL



# unit test (end-to-end)
ifeq ($(MAKECMDGOALS),test)
LDLIBS := $(filter-out -lfftw3f,$(LDLIBS))
endif
test: bin/plambda bin/imprintf
	bin/plambda zero:512x512 "randg randg randl randg randg 5 njoin" \
	| bin/plambda - "split rot del hypot" | bin/imprintf "%s%e%r%i%a\n" \
	| grep -q 3775241.440161.730120.0037888911.8794


# hacks (compatibility flags for old compilers)

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) $(CFLAGS) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif
