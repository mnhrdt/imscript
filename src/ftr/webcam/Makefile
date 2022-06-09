CFLAGS = -g
CFLAGS = -O3 -march=native

LDLIBS = -lX11 -ljpeg

BIN = camview
OBJ = camview.o cam.o ftr.o iio.o

$(BIN) : $(OBJ)

iio.o : CFLAGS += -DI_CAN_HAS_LIBJPEG

clean: ; $(RM) $(BIN) $(OBJ)
