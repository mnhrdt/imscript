# user-editable configuration

# swap these two lines to change between debug and release
CFLAGS = -g -Wall -Wextra
CFLAGS = -march=native -O3 -ffast-math

#TODO: add make options to switch between the X11 and freeglut interfaces

# -- end of user-editable part


# variables
OBJ = iio.o fancy_image.o ftr.o egm96.o
LIB = libiio.a
BIN = viho fpan fpantiff rpcflip icrop

# default LDLIBS contains everything
LDLIBS = $(LIB) -ljpeg -ltiff -lpng -lz -lm -lX11


# rules

# default rule (builds everything)
default: $(BIN)

# rule to build all the individual executables
$(BIN) : $(LIB)
#$(BIN) : % : %.c $(LIB)
#	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

# rule to build the "iio" library interfaces, linkable with -liio
$(LIB): $(OBJ)
	$(AR) rc $@ $^

# automatic generation of dependencies
.depend.mk: $(MAKEFILE_LIST)
	-$(CC) -MM $(BIN:=.c) | sed 's/\.o:/:/' > .depend.mk
-include .depend.mk


# bureaucracy
clean:
	$(RM) $(LIB) $(OBJ) $(BIN) im
depend: .depend.mk
.PHONY: default clean depend
