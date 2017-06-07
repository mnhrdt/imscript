SRC	= $(addsuffix .c, \
	setdim fftshift imdim sum qauto fft sort closeup gblur \
	plambda synflow flowdiv viewflow faxpb downsa component \
	qeasy vecstack backflow rgbcube iminfo fnorm faxpby)
OBJ	= $(SRC:.c=.o)
BIN	= $(SRC:.c=)

COPT	= -std=c99 -pedantic-errors -Wall -Wextra -Wshadow -Wno-unused
CFLAGS	= $(COPT) -g -I../iio/
DESTDIR	= bin
LDFLAGS	= -ltiff -lpng -lfftw3f
LDLIBS	= ../iio/iio.o

default	:	$(BIN)

# cancel default rules
.SUFFIXES	: 
# object compilation
%.o	: %.c
	$(CC) -c $(CFLAGS) $< -o $@
# final link
%	: %.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY	: install clean distclean
install	: $(BIN)
	install $^ $(DESTDIR)/
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN)
