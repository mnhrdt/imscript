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

# default target
default: $(addprefix bin/,$(BIN))

# rules
bin/% : src/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(LIB) : $(LIB)($(OBJ))

# bureaucracy
clean:
	$(RM) $(OBJ) $(LIB) $(BIN)
.PHONY: default clean

# hacks (compatibility hacks for old compilers)
-include src/make.hacks.mk
