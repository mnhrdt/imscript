CFLAGS = -g
CFLAGS = -O3 -march=native

LDLIBS = -lm -lfftw3f -lX11 -ljpeg

BIN = camview corrview
OBJ = cam.o ftr.o iio.o

all : $(BIN)

$(BIN) : $(OBJ)

iio.o : CFLAGS += -DI_CAN_HAS_LIBJPEG

clean: ; $(RM) $(BIN) $(OBJ)
