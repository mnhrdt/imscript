CFLAGS ?= -march=native -O3 -DFTR_BACKEND='f'

LDLIBS = -ljpeg -ltiff -lpng -lz -lm -lGL -lglut
LDLIBS = -ljpeg -ltiff -lpng -lfftw3f -lz -lm -lX11

BIN = viho fpan fpantiff rpcflip icrop powerkill dosdo epiview vnav
OBJ = iio.o fancy_image.o ftr.o egm96.o
LIB = libiio.a

default : $(BIN)
$(BIN)  : $(LIB)
$(LIB)  : $(LIB)($(OBJ))

clean:
	$(RM) $(LIB) $(OBJ) $(BIN)
.PHONY: default clean
