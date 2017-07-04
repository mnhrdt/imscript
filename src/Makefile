# user-editable configuration
CFLAGS ?= -march=native -O3

# required libraries
LDLIBS  = -lfftw3f -lm -ljpeg -ltiff -lpng -lz

# variables
OBJ = iio.o fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac      \
      srmatch tiffu siftu crop lrcat tbcat fftshift imflip bmms registration \
      blur fft dct dht flambda fancy_crop fancy_downsa iion iion_u16

# targets
default    : $(BIN)

# dependencies (manual)
$(BIN) : % : %.o $(OBJ)

# bureaucracy
clean      : ; $(RM) *.a *.o $(BIN)
.PHONY     : default clean

# automatic generation of missing dependencies
deps.mk    : ; $(CC) -MM $(BIN:=.c) > deps.mk
-include deps.mk

# hack for older compilers (adds gnu99 option to CFLAGS, if necessary)
ifeq (,$(shell $(CC) $(CFLAGS) -dM -E - < /dev/null | grep __STDC_VERSION__))
CFLAGS := $(CFLAGS) -std=gnu99
endif
