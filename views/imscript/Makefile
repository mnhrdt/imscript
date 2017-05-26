# user-editable configuration (with reasonable defaults)
CFLAGS ?= -march=native -O3

# required libraries
LDLIBS = -ljpeg -ltiff -lpng -lz -lfftw3f -lm
LDLIBS_FTR = $(LDLIBS) -lX11

# variables
OBJ = src/iio.o src/fancy_image.o
LIB = src/libiio.a
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa iion                       \

OBJ_FTR = $(OBJ) src_ftr/ftr.o src_ftr/egm96.o
LIB_FTR = src_ftr/libftr.a
BIN_FTR = viho fpan fpantiff rpcflip icrop

OBJ_ALL = $(OBJ) $(OBJ_FTR)
LIB_ALL = $(LIB) $(LIB_FTR)
BIN_ALL = $(addprefix bin/,$(BIN) $(BIN_FTR))

# default target
default: $(BIN_ALL)

# rules
bin/% : src/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/% : src_ftr/%.c $(LIB_FTR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS_FTR)

$(LIB)     : $(LIB)($(OBJ))
$(LIB_FTR) : $(LIB_FTR)($(OBJ_FTR))

# bureaucracy
clean:
	$(RM) $(OBJ_ALL) $(LIB_ALL) $(BIN_ALL)
.PHONY: default clean

# hacks
include src/.make.hacks.mk
