# user-editable configuration

# swap these two lines to change between debug and release
CFLAGS = -g -Wall -Wextra
CFLAGS = -march=native -Os

# set the following variables to "yes" to disable some dependences
DISABLE_IMAGE_LIBS = no
DISABLE_FFTW3 = no

# -- end of user-editable part


# variables
OBJ = iio.o fancy_image.o
BIN = plambda vecov morsi downsa upsa ntiply censust dither \
      qauto qeasy \
      homwarp synflow backflow flowinv \
      nnint bdint amle simpois \
      ghisto contihist \
      fontu imprintf pview viewflow flowarrows palette \
      ransac srmatch \
      tiffu siftu \
      crop lrcat tbcat fftshift \
      bmms registration \
      blur fft dct dht \
      flambda fancy_crop fancy_downsa \
      iion

# default LDLIBS contains with everything
IMAGE_LIBS = -ljpeg -ltiff -lpng -lz
LDLIBS = libiio.a $(IMAGE_LIBS) -lfftw3f -lm

# if requested, remove image libraries
ifeq ($(DISABLE_IMAGE_LIBS),yes)
FANCY_IMAGE_FLAGS = -DFANCY_IMAGE_DISABLE_TIFF
IIO_FLAGS         = -DIIO_DISABLE_IMGLIBS
LDLIBS            := $(filter-out $(IMAGE_LIBS),$(LDLIBS))
BIN               := $(filter-out tiffu,$(BIN))
endif

# if requestesd, remove FFT-depending programs
ifeq ($(DISABLE_FFTW3),yes)
LDLIBS := $(filter-out -lfftw3f,$(LDLIBS))
BIN    := $(filter-out blur fft dct dht,$(BIN))
endif

# augment flags with target-specific variables
fancy_image.o : override CFLAGS += $(FANCY_IMAGE_FLAGS)
iio.o         : override CFLAGS += $(IIO_FLAGS)



# rules

# default rule (builds everything)
default: im $(BIN)

# rule to build all the individual executables
$(BIN) : % : %.c libiio.a
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

# rule to build the "iio" library interface, linkable with -liio
libiio.a: $(OBJ)
	$(AR) rc $@ $^

# rule to build a single executable "im" containing all the programs
im: im.c $(BIN:=.c) libiio.a
	@true > all_mains.inc
	@for i in $(BIN); do printf "\tMAIN($$i)\n" >> all_mains.inc ; done
	$(CC) $(CFLAGS) -DHIDE_ALL_MAINS im.c $(BIN:=.c) $(LDLIBS) -o $@


# automatic generation of dependencies
.depend.mk: $(MAKEFILE_LIST)
	-$(CC) -MM $(BIN:=.c) | sed 's/\.o:/:/' > .depend.mk
-include .depend.mk



# bureaucracy
clean:
	$(RM) *.a *.o $(BIN) im
list_bin:
	@echo $(BIN)
depend: .depend.mk
.PHONY: default clean depend list_bin
