# user-editable configuration
CFLAGS ?= -march=native -Os

# variables
OBJ = iio.o fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac      \
      srmatch tiffu siftu crop lrcat tbcat fftshift imflip bmms registration \
      blur fft dct dht flambda fancy_crop fancy_downsa iion

# libraries
LDLIBS = -lfftw3f -lm -ljpeg -ltiff -lpng -lz

# default target
default: $(BIN)

# rule to build each program
$(BIN) : % : %.o $(OBJ)

# automatic generation of dependences
depend.mk:
	$(CC) -MM $(BIN:=.c) > depend.mk
include depend.mk

# bureaucracy
clean:
	$(RM) *.a *.o $(BIN)
.PHONY: default clean
