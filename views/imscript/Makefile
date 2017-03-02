SHELL = bash -O extglob

export CFLAGS = -march=native -O3

all: src
	$(MAKE) -C src -j 3
	$(MAKE) -C src_ftr -j 3
	cp src/+([a-z]) bin
	cp src_ftr/+([a-z]) bin

clean:
	$(MAKE) -C src clean
	$(MAKE) -C src_ftr clean
	$(RM) bin/*
