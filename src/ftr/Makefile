CFLAGS ?= -march=native -O3
#CFLAGS ?= -g -fsanitize=address

LDLIBS  = -ljpeg -ltiff -lpng -lfftw3f -lz -lm -lX11

BIN = viho fpan fpantiff rpcflip icrop powerkill dosdo epiview vnav \
      pleview srt fpanflip wifpan
OBJ = iio.o ftr.o egm96.o

default : $(BIN)
$(BIN)  : $(OBJ)

clean   : ; $(RM) $(OBJ) $(BIN)
.PHONY  : default clean
