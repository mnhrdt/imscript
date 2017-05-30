CFLAGS ?= -march=native -Os
OBJ = iio.o fancy_image.o
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac      \
      srmatch tiffu siftu crop lrcat tbcat fftshift imflip bmms registration \
      blur fft dct dht flambda fancy_crop fancy_downsa iion
LDLIBS = -lfftw3f -lm -ljpeg -ltiff -lpng -lz
default: $(BIN)
$(BIN) : % : %.o $(OBJ)
depend.mk:
	$(CC) -MM $(BIN:=.c) > depend.mk
clean:
	$(RM) *.a *.o $(BIN)
.PHONY: clean default
-include depend.mk
