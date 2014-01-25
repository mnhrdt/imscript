CC = c99

LFLAGS = -ltiff -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -g
CFLAGS = -O3 -Wall -Wextra -pedantic

all: cut3 join3 translation quantize registration

cut3: cut3.c iio.c
	$(CC) $(CFLAGS) iio.c cut3.c $(LFLAGS) -o cut3

join3: join3.c iio.c
	$(CC) $(CFLAGS) iio.c join3.c $(LFLAGS) -o join3

translation: translation.c iio.c
	$(CC) $(CFLAGS) iio.c translation.c $(LFLAGS) -o translation

registration: registration.c iio.c
	$(CC) $(CFLAGS) iio.c registration.c $(LFLAGS) -o registration

quantize: quantize.c iio.c
	$(CC) $(CFLAGS) iio.c quantize.c $(LFLAGS) -o quantize

clean:
	rm -f cut3 join3 translation quantize registration
