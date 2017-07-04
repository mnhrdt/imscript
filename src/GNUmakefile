CFLAGS ?= -march=native -O3
LDLIBS  = -lfftw3f -lm -ljpeg -ltiff -lpng -lz
OBJ = iio.o fancy_image.o
BIN = $(shell cat TARGETS)
default    : $(BIN)
$(BIN)     : $(OBJ)
clean      : ; $(RM) *.o $(BIN)
deps.mk    : ; $(CC) -MM $(BIN:=.c) > deps.mk
-include deps.mk
