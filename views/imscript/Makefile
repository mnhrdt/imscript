export CFLAGS = -march=native -O3

# TODO: put some order...
#
# there are a few types of files here
#
# 1. header-only files like "abstract_heap.h" or "smapa.h"
# These are never compiled alone and pose no problem for the build system.
# They only need to be at a place where they are found by #include
#
# 2. C files that contain a "main_*" function:
#
# 2.1. depending on no libraries whatsoever : LDLIBS=""
# 2.2. depending on iio.o and associated stuff : LDLIBS+="-lpng..."
# 2.3. like 2.2., but depends also on fancy_image.o
# 2.4. like 2.2., but depends also on fftw3 : LDLIBS+="-lfftw3f"
# 2.5. like 2.2., but depends also on gsl : LDLIBS+="-lgsl -lgslcblas"
# 2.6. interactive programs that depend on ftr.o and X11 or freeglut
#
# All the files of type "2" can be either compiled standalone, or #included
# inside another file.
#
# 3. C files that do not contain a main_ function
#
# The files of type "3" are intended to be included into other programs, and
# all their functions are marked static.
#
# 4.  Only three files are intended to be compiled as .o
#
# 	4.1. iio.o
# 	4.2. fancy_image.o
# 	4.3. ftr.o
#
# For other files, this is impossible: if you want to use their functions, just
# include the ".c".  The traditional ".h/.c" duality is most of the time an
# ugly overkill.
#
#

all:
	$(MAKE) -C src -j 3
	@for i in `$(MAKE) -s -C src list_bin`; do cp src/$$i bin ; done

clean:
	$(MAKE) -C src clean
	$(RM) bin/*
