CFLAGS ?= -march=native -O3                     # user-editable config
LDLIBS  = -lfftw3f -lm -ljpeg -ltiff -lpng -lz  # required libraries
OBJ     = iio.o fancy_image.o                   # object files
BIN     = $(shell cat TARGETS)                  # executables
default : $(BIN)                                # default target
$(BIN)  : $(OBJ)                                # dependencies
clean   : ; $(RM) *.o $(BIN)                    # bureaucracy
deps.mk : ; $(CC) -MM $(BIN:=.c) > deps.mk      # retrieve missing dependencies
-include deps.mk                                # include missing dependencies
