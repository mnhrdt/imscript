# user-editable configuration
CFLAGS ?= -march=native -Os

# required libraries
LDLIBS = -ljpeg -ltiff -lpng -lz -lm

# targets
OBJ = iio.o fancy_image.o
LIB = libiio.a
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview palette viewflow flowarrows ransac blur \
      srmatch tiffu siftu crop fftshift imflip lrcat tbcat bmms registration \
      fft dct dht flambda fancy_crop fancy_downsa iion                       \

# dependences
default : $(BIN)
$(BIN)  : $(LIB)
$(LIB)  : $(LIB)($(OBJ))

# bureaucracy
clean:
	$(RM) $(LIB) $(OBJ) $(BIN)
.PHONY: default clean
