CC = c99

LFLAGS = -ltiff -lpng -lm
WFLAGS = -Wall -Wextra

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = $(WFLAGS) -g
CFLAGS = $(WFLAGS) -Ofast

PROGRAMS = cut3 join3 translation quantize registration tiles paircrop combine_homographies viewflow

all: $(PROGRAMS)

iio.o: iio.c
	$(CC) $(CFLAGS) -c iio.c

%: %.c iio.o
	$(CC) $(CFLAGS) $^ $(LFLAGS) -o $@

clean:
	rm -f $(PROGRAMS) iio.o
