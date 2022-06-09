CFLAGS = -g
CFLAGS = -O3 -march=native -DI_CAN_HAS_LIBJPEG

camview: camview.c cam.c cam.h ftr.c iio.c
	$(CC) $(CFLAGS) camview.c cam.c ftr.c iio.c -lX11 -ljpeg -o $@

clean: ; $(RM) camview
