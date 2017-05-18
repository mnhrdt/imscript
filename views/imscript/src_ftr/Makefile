CFLAGS ?= -march=native -O3 -DFTR_BACKEND='f'

#LDLIBS = -ljpeg -ltiff -lpng -lz -lm -lX11
LDLIBS = -ljpeg -ltiff -lpng -lz -lm -lGL -lglut

BIN = viho fpan fpantiff rpcflip icrop
OBJ = iio.o fancy_image.o ftr.o egm96.o
LIB = libiio.a

default : $(BIN)
$(BIN)  : $(LIB)
$(LIB)  : $(LIB)($(OBJ))

clean:
	$(RM) $(LIB) $(OBJ) $(BIN)
.PHONY: default clean
