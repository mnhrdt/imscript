CFLAGS ?= -march=native -O3
#CFLAGS ?= -g -fsanitize=address

LDLIBS  = -ljpeg -ltiff -lpng -lfftw3f -lz -lm -lX11

BIN = $(shell cat TARGETS)
OBJ = iio.o fancy_image.o ftr.o egm96.o

default : $(BIN)
$(BIN)  : $(OBJ)

clean   : ; $(RM) $(OBJ) $(BIN)
.PHONY  : default clean
