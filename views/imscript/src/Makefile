# setup
LDLIBS = -ljpeg -ltiff -lpng -lz -lm
LDLIBS = -ljpeg -ltiff -lpng -lfftw3f -lz -lm
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
      iion

# rules
default: $(BIN) #im

$(BIN) : % : %.c iio.o            # individual executables
	$(CC) $(CFLAGS) $< iio.o $(LDLIBS) -o $@

im: im.c $(BIN:=.c) iio.o         # one static box with all the stuff inside
	@true > all_mains.inc
	@for i in $(BIN); do printf "\tMAIN($$i);\n" >> all_mains.inc ; done
	$(CC) $(CFLAGS) -DHIDE_ALL_MAINS $^ $(LDLIBS) -o $@

# automatic generation of dependencies
.depend.mk: $(MAKEFILE_LIST)
	-$(CC) -MM $(BIN:=.c) | sed 's/\.o:/:/' > .depend.mk
-include .depend.mk

# bureaucracy
clean:
	$(RM) *.o $(BIN) im
list_bin:
	@echo $(BIN)
depend: .depend.mk
.PHONY: default clean depend list_bin
