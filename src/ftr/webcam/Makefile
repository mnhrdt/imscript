CFLAGS = -g
CFLAGS = -O3 -march=native

camview: camview.c cam.c cam.h ftr.c
	$(CC) $(CFLAGS) camview.c cam.c ftr.c -lX11 -o $@

clean: ; $(RM) camview
