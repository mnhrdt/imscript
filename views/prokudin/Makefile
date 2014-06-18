CC = c99

LFLAGS = -ltiff -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -g
CFLAGS = -O3 -Wall -Wextra

PROGRAMS = cut3 join3 translation quantize registration

all: $(PROGRAMS)

iio.o: iio.c
	$(CC) $(CFLAGS) -c iio.c

%: %.c iio.o
	$(CC) $(CFLAGS) $^ $(LFLAGS) -o $@

clean:
	rm -f $(PROGRAMS) iio.o
