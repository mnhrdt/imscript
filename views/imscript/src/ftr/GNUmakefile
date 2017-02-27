# configuration
CFLAGS = -g -Wall -Wextra
CFLAGS = -march=native -O3 -ffast-math

# variables
OBJ = iio.o fancy_image.o ftr.o egm96.o
LIB = libiio.a
BIN = viho fpan fpantiff rpcflip icrop
LDLIBS = -ljpeg -ltiff -lpng -lz -lm -lX11

# dependencies
default: $(BIN)
$(BIN) : $(LIB)
$(LIB) : $(LIB)($(OBJ))

# bureaucracy
clean:
	$(RM) $(LIB) $(OBJ) $(BIN)
.PHONY: default clean
