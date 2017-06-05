CFLAGS ?= -march=native -Os
LDLIBS  = -lfftw3f -lgsl -lgslcblas -lm -ljpeg -ltiff -lpng -lz
OBJ = iio.o fancy_image.o
BIN := $(shell cat all_mains)
default    : $(BIN)
$(BIN) : % : %.o $(OBJ)
deps.mk    : ; $(CC) -MM $(BIN:=.c) > deps.mk
clean      : ; $(RM) *.a *.o $(BIN)
-include deps.mk
