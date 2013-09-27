CC = c99

LFLAGS = -ltiff -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -g
CFLAGS = -O3 -fopenmp

SRC = centroid.c iio.c optical_flow/pyramidal_horn_schunck.c optical_flow/smooth_inpainting.c optical_flow/optical_flow.c

all: centroid combine

centroid: $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LFLAGS) -o centroid

combine: combine.c iio.c
	$(CC) $(CFLAGS) combine.c iio.c $(LFLAGS) -o combine

clean:
	rm -f centroid combine
