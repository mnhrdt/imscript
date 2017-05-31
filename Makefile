CFLAGS ?= -march=native -O3
LDLIBS = -ljpeg -ltiff -lpng -lz -lfftw3f -lm

OBJ = src/iio.o src/fancy_image.o
LIB = src/libiio.a
BIN = plambda vecov veco vecoh morsi downsa upsa ntiply censust dither qauto \
      qeasy homwarp synflow backflow flowinv nnint bdint amle simpois ghisto \
      contihist fontu imprintf pview viewflow flowarrows palette ransac blur \
      srmatch tiffu siftu crop lrcat tbcat fftshift bmms registration imflip \
      fft dct dht flambda fancy_crop fancy_downsa iion                       \

BIN := $(addprefix bin/,$(BIN))

default: $(BIN)

bin/%  : src/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(LIB) : $(LIB)($(OBJ))




# END OF CORE PART, THE REST OF THIS FILE IS NOT ESSENTIAL

# unit test (end-to-end)
ifeq ($(MAKECMDGOALS),test)
LDLIBS := $(filter-out -lfftw3f,$(LDLIBS))
endif
test: bin/plambda bin/imprintf
	bin/plambda zero:512x512 "randg randg randl randg randg 5 njoin" \
	| bin/plambda - "split rot del hypot" | bin/imprintf "%s%e%r%i%a\n" \
	| grep -q 3775241.440161.730120.0037888911.8794


# hack (compatibility flags for old compilers)

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) $(CFLAGS) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif




# exotic targets, not built by default
# FTR: interactive tools, require X11 or freeglut
# LEG: "legacy" tools, or those requiring GSL

OBJ_FTR = $(OBJ) src/ftr/ftr.o src/ftr/egm96.o
LIB_FTR = src/ftr/libftr.a

BIN_FTR = viho fpan fpantiff rpcflip icrop powerkill dosdo epiview vnav
BIN_LEG = $(shell cat src/legacy/all_mains)

BIN_FTR := $(addprefix bin/,$(BIN_FTR))
BIN_LEG := $(addprefix bin/,$(BIN_LEG))

LDLIBS_FTR = $(LDLIBS) -lX11
LDLIBS_LEG = $(LDLIBS) -lgsl -lgslcblas

OBJ_ALL = $(OBJ) $(OBJ_FTR)
LIB_ALL = $(LIB) $(LIB_FTR)
BIN_ALL = $(BIN) $(BIN_FTR) $(BIN_LEG)

full  : default ftr legacy
ftr   : $(BIN_FTR)
legacy: $(BIN_LEG)


bin/% : src/ftr/%.c $(LIB_FTR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS_FTR)

bin/% : src/legacy/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS_LEG)

$(LIB_FTR) : $(LIB_FTR)($(OBJ_FTR))



# bureaucracy
clean:
	$(RM) $(OBJ_ALL) $(LIB_ALL) $(BIN_ALL)
.PHONY: default full ftr legacy clean
