CFLAGS ?= -march=native -Os
LDLIBS  = -lfftw3f -lm -ljpeg -ltiff -lpng -lz
OBJ = iio.o fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac      \
      srmatch tiffu siftu crop lrcat tbcat fftshift imflip bmms registration \
      blur fft dct dht flambda fancy_crop fancy_downsa iion
default    : $(BIN)
$(BIN) : % : %.o $(OBJ)
deps.mk    : ; $(CC) -MM $(BIN:=.c) > deps.mk
clean      : ; $(RM) *.a *.o $(BIN)
-include deps.mk
