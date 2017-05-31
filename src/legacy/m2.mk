CFLAGS ?= -march=native -Os
LDLIBS  = -lfftw3f -lm -ljpeg -ltiff -lpng -lz
OBJ = iio.o fancy_image.o
BIN := $(shell for i in *.c;do cc -E $$i|grep -q '^int main.int '&&echo $$i;done|sed 's/\.c//')
default    : $(BIN)
$(BIN) : % : %.o $(OBJ)
deps.mk    : ; $(CC) -MM $(BIN:=.c) > deps.mk
clean      : ; $(RM) *.a *.o $(BIN)
-include deps.mk
