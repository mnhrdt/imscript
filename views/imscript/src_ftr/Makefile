CFLAGS ?= -march=native -O3

LDLIBS = -ljpeg -ltiff -lpng -lz -lm -lX11

BIN = viho fpan fpantiff rpcflip icrop
OBJ = iio.o fancy_image.o ftr.o egm96.o
LIB = libiio.a

default : $(BIN)
$(BIN)  : $(LIB)
$(LIB)  : $(LIB)($(OBJ))

clean:
	$(RM) $(LIB) $(OBJ) $(BIN)
.PHONY: default clean
