# user-editable configuration with reasonable defaults
CFLAGS ?= -march=native -Os

# libraries
LDLIBS = -ljpeg -ltiff -lpng -lz -lfftw3f -lm

# targets
OBJ = iio.o fancy_image.o
LIB = libiio.a
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa iion                       \

# dependences
default: $(BIN)
$(BIN) : $(LIB)
$(LIB) : $(LIB)($(OBJ))

# bureaucracy
clean:
	$(RM) $(OBJ) $(LIB) $(BIN)
.PHONY: default clean

# hacks (optional, needed only to solve bizarre problems on old systems)
include .make.hacks.mk

