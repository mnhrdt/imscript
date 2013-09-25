CC = c99

LFLAGS = -ltiff -lpng -lm -fopenmp

CFLAGS = -O3
CFLAGS = -g

SRC = centroid.c iio.c optical_flow/pyramidal_horn_schunck.c optical_flow/smooth_inpainting.c optical_flow/optical_flow.c

centroid: $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LFLAGS) -o centroid

clean:
	rm -f centroid
